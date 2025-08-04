use self::{
    av_oracle::AVOracle,
    enums::{AVError, OracleState},
};

use icicle_core::curve::{Curve,Affine};

mod av_oracle;
pub(crate) mod enums;
/////////////////////////////////////////////////////////////////////////////
// This part directly uses the original code to ensure compatibility.
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

    pub fn register_msgs(&mut self, msgs: &Vec<Affine::<C>>, pos: usize) -> Result<(), AVError> {
        assert_eq!(msgs.len(), N);
        match self.state {
            OracleState::Round1Ongoing => {
                for i in 0..N {
                    self.av_oracles[i].register_msg(msgs[i], pos)?;
                }

                self.first_msgs_registered += 1;
                if self.first_msgs_registered == B {
                    self.state = OracleState::Round1Completed;
                }
            }
            OracleState::Round1Completed => {
                println!("in here");
                return Err(AVError::WrongState(String::from(
                    "Auctioneer: Message can't be sent in Round1Completed state",
                )));
            }
            OracleState::Round2Ongoing => {
                for i in 0..N {
                    self.av_oracles[i].register_msg(msgs[i], pos)?;
                }

                self.second_msgs_registered += 1;
                if self.second_msgs_registered == B {
                    self.state = OracleState::Round2Completed;
                }
            }
            OracleState::Round2Completed => {
                return Err(AVError::WrongState(String::from(
                    "Auctioneer: Message can't be sent in Round2Completed state",
                )))
            }
            OracleState::Completed => {
                return Err(AVError::WrongState(String::from(
                    "Auctioneer: Message can't be sent in Completed state",
                )))
            }
        }

        Ok(())
    }
////////////
/////////////////////////////////////////////////////////////////////////////
    pub fn output_first_round(&mut self) -> Vec<Vec<Affine::<C>>> {
        assert_eq!(self.state, OracleState::Round1Completed);
        self.state = OracleState::Round2Ongoing;
        
        // may influence efficiency
        let mut result = Vec::new();
        for av_i in &mut self.av_oracles {
            result.push(av_i.output_first_round());
        }
        result
    }

    pub fn output_second_round(&mut self) -> Vec<Affine::<C>> {
        assert_eq!(self.state, OracleState::Round2Completed);
        self.state = OracleState::Completed;

        // may influence efficiency
        let mut result = Vec::new();
        for av_i in &mut self.av_oracles {
            result.push(av_i.output_second_round());
        }
        result
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

        let mut a = Auctioneer::<N, B, Bn254CurveCfg>::new();
        let mut secrets = vec![vec![<Bn254CurveCfg as Curve>::ScalarField::zero(); N]; B];
        let mut first_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];

        // initialize n msgs for each party
        for i in 0..B {
            for j in 0..N {
                secrets[i][j] = ScalarCfg::generate_random(1)[0];
            }
        }

        // initialize n msgs for each party
        for i in 0..B {
            for j in 0..N {
                first_msgs[i][j] = g.mul(secrets[i][j]).into();
            }
        }

        // each party sends it's first round msgs
        for i in 0..B {
            a.register_msgs(&first_msgs[i], i).unwrap();
        }

        // we get output for each party per round
        // where each row is of len B (output of av for each party)
        let fr_result = a.output_first_round();

        let mut second_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];
        for i in 0..B {
            for j in 0..N {
                second_msgs[i][j] = Bn254CurveCfg::mul_scalar(fr_result[j][i].into(), secrets[i][j]).into();
            }
        }

        // each party sends it's second round msgs
        for i in 0..B {
            a.register_msgs(&second_msgs[i], i).unwrap();
        }

        let av_results = a.output_second_round();
        // since each party used same secret, expect that all outputs are 0
        for i in 0..N {
            assert_eq!(av_results[i], Affine::<Bn254CurveCfg>::zero());
        }
    }
}