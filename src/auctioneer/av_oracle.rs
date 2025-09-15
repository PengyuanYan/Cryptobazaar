//! Anonymous Veto per‑position oracle.
//!
//! An `AVOracle` encapsulates the state machine and arithmetic for a single
//! auction position/bucket. `Auctioneer` push bidder messages
//! into the oracle in two phases, and then extract the round outputs.
use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::traits::FieldImpl;
use super::enums::{AVError, OracleState};
use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;
use crate::utils::{msm_gpu, my_msm, get_coeffs_of_poly, get_device_is_cpu_or_gpu, to_affine_batched};

// this structure is directly used from the original code
#[derive(Clone)]
pub struct AVOracle<const B: usize, C: Curve> {
    state: OracleState,
    first_msgs: Vec<Affine::<C>>,
    first_msgs_registered: usize,
    second_msgs: Vec<Affine::<C>>,
    second_msgs_registered: usize,
}

impl<const B: usize, C: Curve + icicle_core::msm::MSM<C>> AVOracle<B, C> {
    pub fn new() -> Self {
        Self {
            state: OracleState::Round1Ongoing,
            first_msgs: vec![Affine::<C>::zero(); B],
            first_msgs_registered: 0,
            second_msgs: vec![Affine::<C>::zero(); B],
            second_msgs_registered: 0,
        }
    }

    /// Returns the message bucket for the current round, or an error if sending is not allowed.
    fn current_round_msgs(&self) -> Result<&[Affine<C>], AVError> {
        let wrong = |s: &str| AVError::WrongState(format!("AV: Message can't be sent in {s} state"));

        match self.state {
            OracleState::Round1Ongoing => Ok(&self.first_msgs),
            OracleState::Round2Ongoing => Ok(&self.second_msgs),
            OracleState::Round1Completed => Err(wrong("Round1Completed")),
            OracleState::Round2Completed => Err(wrong("Round2Completed")),
            OracleState::Completed       => Err(wrong("Completed")),
        }
    }

    fn msg_validity(&self, msg: Affine<C>, pos: usize) -> Result<(), AVError> {
        if pos >= B {
            return Err(AVError::WrongPosition(format!("AV: Position must be < {}", B)));
        }

        if msg == Affine::<C>::zero() {
            return Err(AVError::WrongMsg(String::from("AV: Message can't be zero")));
        }

        // Ensure this position isn't already used for the active round.
        let msgs = self.current_round_msgs()?;
        if msgs[pos] != Affine::<C>::zero() {
            return Err(AVError::MsgAlreadySent(pos));
        }

        Ok(())
    }
    
    /// Register an incoming bidder message, validating the current phase.
    pub fn register_msg(&mut self, msg: Affine::<C>, pos: usize) -> Result<(), AVError> {
        self.msg_validity(msg, pos)?;
        // Route by phase and reject messages that don’t belong to the current phase.
        match self.state {
            OracleState::Round1Ongoing => {
                self.first_msgs[pos] = msg;
                self.first_msgs_registered += 1;

                if self.first_msgs_registered == B {
                    self.state = OracleState::Round1Completed;
                }
            }
            OracleState::Round2Ongoing => {
                self.second_msgs[pos] = msg;
                self.second_msgs_registered += 1;

                if self.second_msgs_registered == B {
                    self.state = OracleState::Round2Completed;
                }
            }
            _ => panic!("AV: Something went wrong in register_msg"),
        }

        Ok(())
    }
    
    /// The core logic of anonymous veto round 1 and return the per‑position first‑round output.
    /// It exploit the trick mentioned in Cryptobazaar.
    pub fn output_first_round(&mut self, cpu_or_gpu: usize) -> Vec<Affine::<C>> {
        assert_eq!(self.state, OracleState::Round1Completed);

        let batched = 0; // if use the to_affine_batched function. 0 for no, 1 for yes.

        let mut x: Vec<C::ScalarField> = (0..B)
            .map(|i| if i == B - 1 { C::ScalarField::zero() }
                     else { C::ScalarField::one() })
        .collect();
        
        if batched == 1 {
            let mut outputs = vec![Affine::<C>::zero(); B];
            let projective_output = my_msm(&x, &self.first_msgs, cpu_or_gpu);
            outputs[B-1] = projective_output.into();

            /*
            b1 0 -1 -1 -1
            b2 1  0 -1 -1
            b3 1  1  0 -1
            b4 1  1  1  0
            */
            for i in 0..(B - 1) {
                let idx = B - 2 - i;
                let projective_output = outputs[idx + 1].to_projective() - self.first_msgs[idx + 1].to_projective() - self.first_msgs[idx].to_projective();
                outputs[idx] = projective_output.into();
            }

            self.state = OracleState::Round2Ongoing;
            
            return outputs;
        } else {
            let mut outputs = vec![Projective::<C>::zero(); B];
            let projective_output = my_msm(&x, &self.first_msgs, cpu_or_gpu);
            outputs[B-1] = projective_output;

            for i in 0..(B - 1) {
                let idx = B - 2 - i;
                let projective_output = outputs[idx + 1] - self.first_msgs[idx + 1].into() - self.first_msgs[idx].into();
                outputs[idx] = projective_output;
            }

            self.state = OracleState::Round2Ongoing;

            let outputs = to_affine_batched(&outputs);
            return outputs;
        }
    }
    
    /// Fininsh round 2 and return the per‑position second‑round result.
    /// This one should used with the parallel caller.
    pub fn output_second_round_parallel(&mut self, cpu_or_gpu: usize) -> Projective::<C> {
        assert_eq!(self.state, OracleState::Round2Completed);
        self.state = OracleState::Completed;

        let ones = vec![C::ScalarField::one(); B];

        my_msm(&ones, &self.second_msgs, cpu_or_gpu)
    }
    
    /// Fininsh round 2 and return the per‑position second‑round result.
    /// This one should used with the not parallel caller.
    pub fn output_second_round_not_parallel(&mut self, cpu_or_gpu: usize) -> Affine::<C> {
        assert_eq!(self.state, OracleState::Round2Completed);
        self.state = OracleState::Completed;

        let ones = vec![C::ScalarField::one(); B];

        my_msm(&ones, &self.second_msgs, cpu_or_gpu).into()
    }
}

#[cfg(test)]
mod av_oracle_tests {
    use std::ops::Mul;

    use super::AVOracle;
    const B: usize = 1024;

    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg};
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::traits::GenerateRandom;

    use icicle_bn254::curve::ScalarCfg;

    #[test]
    fn no_veto() {
        let g = Bn254CurveCfg::get_generator();

        let mut av_oracle = AVOracle::<B, Bn254CurveCfg>::new();

        let party_secrets: Vec<_> = (0..B).map(|_| ScalarCfg::generate_random(1)[0]).collect();
        
        let mut first_msgs = Vec::with_capacity(party_secrets.len());

        for mi in &party_secrets {
            let point = g.mul(*mi);
            first_msgs.push(point.into());
        }

        for i in 0..B {
            av_oracle.register_msg(first_msgs[i], i).unwrap();
        }

        let round_outputs = av_oracle.output_first_round(0);
        let second_msgs: Vec<Affine::<Bn254CurveCfg>> = party_secrets
            .iter()
            .zip(round_outputs.iter())
            .map(|(mi, &ri)| Bn254CurveCfg::mul_scalar(ri.into(), *mi).into())
            .collect();

        for i in 0..B {
            av_oracle.register_msg(second_msgs[i], i).unwrap();
        }

        let output = av_oracle.output_second_round_not_parallel(0);
        assert_eq!(output, Affine::<Bn254CurveCfg>::zero().into());
    }

    #[test]
    fn veto() {
        let g = Bn254CurveCfg::get_generator();

        let mut av_oracle = AVOracle::<B, Bn254CurveCfg>::new();

        let mut party_secrets: Vec<_> = (0..B).map(|_| ScalarCfg::generate_random(1)[0]).collect();

        let mut first_msgs = Vec::with_capacity(party_secrets.len());

        for mi in &party_secrets {
            let point = g.mul(*mi);
            first_msgs.push(point.into());
        }

        for i in 0..B {
            av_oracle.register_msg(first_msgs[i], i).unwrap();
        }

        let round_outputs = av_oracle.output_first_round(0);

        // the case that there is at least one party picks a new secret
        party_secrets[0] = ScalarCfg::generate_random(1)[0];

        let second_msgs: Vec<Affine::<Bn254CurveCfg>> = party_secrets
            .iter()
            .zip(round_outputs.iter())
            .map(|(mi, &ri)| Bn254CurveCfg::mul_scalar(ri.into(), *mi).into())
            .collect();

        for i in 0..B {
            av_oracle.register_msg(second_msgs[i], i).unwrap();
        }

        let output = av_oracle.output_second_round_not_parallel(0);
        assert_ne!(output, Affine::<Bn254CurveCfg>::zero().into());
    }
}