//! Auctioneer is the driver for the anonymous-veto (AV) protocol.
//!
//! This module coordinates `N` independent AV instances across `B` bidders
//! over an elliptic curve `C`. It drives both rounds of the AV protocol,
//! validates state transitions, fans messages out to the per-position oracles,
//! and gathers round outputs.
//!
//! # States
//! The finite-state machine is tracked via [`OracleState`]:
//! - `Round1Ongoing`  -> accepting first-round messages
//! - `Round1Completed`-> first-round outputs available
//! - `Round2Ongoing`  -> accepting second-round messages
//! - `Round2Completed`-> outputs computed; caller can read them
//! - `Completed`      -> terminal state
//! The states are fromt the origainal code.
use self::{
    av_oracle::AVOracle,
    enums::{AVError, OracleState},
};

use crate::utils::get_device_is_cpu_or_gpu;
use icicle_core::curve::{Curve,Affine,Projective};

use rayon::prelude::*;
use crate::utils::to_affine_batched;

mod av_oracle;
pub(crate) mod enums;

#[derive(Clone)]
pub struct Auctioneer<const N: usize, const B: usize, C: Curve + icicle_core::msm::MSM<C>> {
    state: OracleState,
    first_msgs_registered: usize,
    second_msgs_registered: usize,
    av_oracles: Vec<AVOracle<B, C>>,
}

impl<const N: usize, const B: usize, C: Curve + icicle_core::msm::MSM<C>> Auctioneer<N, B, C> {
    pub fn new() -> Self {
        Self {
            state: OracleState::Round1Ongoing,
            first_msgs_registered: 0,
            second_msgs_registered: 0,
            av_oracles: vec![AVOracle::new(); N],
        }
    }
    
    // Small helper that forwards the caller's `msgs`
    // to each per-position oracle, tagging them with `pos`.
    pub fn register_msgs(&mut self, msgs: &[Affine<C>], pos: usize) -> Result<(), AVError> {
        assert_eq!(msgs.len(), N, "expected {N} messages, got {}", msgs.len());

        #[inline]
        fn register_into_oracles<const N: usize, const B: usize, C: Curve + icicle_core::msm::MSM<C>>(
            oracles: &mut [AVOracle<B, C>],
            msgs: &[Affine<C>],
            pos: usize,
        ) -> Result<(), AVError> {
            for (i, msg) in msgs.iter().copied().enumerate() {
                oracles[i].register_msg(msg, pos)?;
            }
            Ok(())
        }

        match self.state {
            OracleState::Round1Ongoing => {
                register_into_oracles::<N, B, C>(&mut self.av_oracles, msgs, pos)?;
                self.first_msgs_registered += 1;
                if self.first_msgs_registered == B {
                    self.state = OracleState::Round1Completed;
                }
                Ok(())
            }
            OracleState::Round2Ongoing => {
                register_into_oracles::<N, B, C>(&mut self.av_oracles, msgs, pos)?;
                self.second_msgs_registered += 1;
                if self.second_msgs_registered == B {
                    self.state = OracleState::Round2Completed;
                }
                Ok(())
            }
            OracleState::Round1Completed => Err(AVError::WrongState(
                "Auctioneer: Message can't be sent in Round1Completed state".into(),
            )),
            OracleState::Round2Completed => Err(AVError::WrongState(
                "Auctioneer: Message can't be sent in Round2Completed state".into(),
            )),
            OracleState::Completed => Err(AVError::WrongState(
                "Auctioneer: Message can't be sent in Completed state".into(),
            )),
        }
    }
    
    /// Finish round 1 and return the vector of first-round outputs per position.
    pub fn output_first_round(&mut self) -> Vec<Vec<Affine::<C>>> {
        assert_eq!(self.state, OracleState::Round1Completed);

        // 1 = run in parallel, 0 = run sequentially
        let parallel = 1; 
        let cpu_or_gpu = get_device_is_cpu_or_gpu();

        self.state = OracleState::Round2Ongoing;

        if parallel == 1 {
        // run in parallel using Rayon
            self.av_oracles
                .par_iter_mut()
                .map(|av| av.output_first_round(cpu_or_gpu))
                .collect()
        } else {
        // run sequentially
            let mut result = Vec::new();
            for av_i in &mut self.av_oracles {
                result.push(av_i.output_first_round(cpu_or_gpu));
            }
            result
        }
    }
    
    /// Finish round 2 and return the aggregated AV outputs.
    pub fn output_second_round(&mut self) -> Vec<Affine::<C>> {
        assert_eq!(self.state, OracleState::Round2Completed);
        self.state = OracleState::Completed;
        
        // 1 = run in parallel, 0 = run sequentially
        let parallel = 1;
        let cpu_or_gpu = get_device_is_cpu_or_gpu();
        
        if parallel == 1 {
        // run in parallel using Rayon
            let projective_outputs = self.av_oracles
                                         .par_iter_mut()
                                         .map(|av| av.output_second_round_parallel(cpu_or_gpu))
                                         .collect::<Vec<Projective<C>>>();
            return to_affine_batched(&projective_outputs)
        } else {
        // run sequentially
            let mut result = Vec::new();
            for av_i in &mut self.av_oracles {
                result.push(av_i.output_second_round_not_parallel(cpu_or_gpu));
            }
            result
        }
    }
}

#[cfg(test)]
mod auctioneer_tests {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_core::curve::{Curve,Affine};
    use icicle_core::traits::FieldImpl;
    use icicle_core::traits::GenerateRandom;
    use std::ops::Mul;
    use icicle_bn254::curve::ScalarCfg;

    use super::Auctioneer;

    const N: usize = 128;
    const B: usize = 32;

    #[test]
    // #[ignore]
    fn test_many_vetos() {
        let g = Bn254CurveCfg::get_generator();

        let mut auctioneer = Auctioneer::<N, B, Bn254CurveCfg>::new();
        let mut secrets = vec![vec![<Bn254CurveCfg as Curve>::ScalarField::zero(); N]; B];
        let mut first_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];

        // generate n random message for each party
        for i in 0..B {
            for j in 0..N {
                secrets[i][j] = ScalarCfg::generate_random(1)[0];
            }
        }

        for i in 0..B {
            for j in 0..N {
                first_msgs[i][j] = g.mul(secrets[i][j]).into();
            }
        }
        
        // send the first round message
        for i in 0..B {
            auctioneer.register_msgs(&first_msgs[i], i).unwrap();
        }
        
        // Outputs of the first av round.
        let fr_result = auctioneer.output_first_round();

        let mut second_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];
        for i in 0..B {
            for j in 0..N {
                second_msgs[i][j] = Bn254CurveCfg::mul_scalar(fr_result[j][i].into(), secrets[i][j]).into();
            }
        }

        // send the second round message
        for i in 0..B {
            auctioneer.register_msgs(&second_msgs[i], i).unwrap();
        }
        
        // all the outputs should be zero
        let av_results = auctioneer.output_second_round();
        for i in 0..N {
            assert_eq!(av_results[i], Affine::<Bn254CurveCfg>::zero());
        }
    }
}