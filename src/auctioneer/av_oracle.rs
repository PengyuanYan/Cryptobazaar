// Oracle that computes outputs of the rounds of AV
use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::traits::FieldImpl;
use super::enums::{AVError, OracleState};
use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;

/////////////////////////////////////////////////////////////////////////////
// This part directly uses the original code to ensure compatibility.
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

    fn msg_validity(&self, msg: Affine::<C>, pos: usize) -> Result<(), AVError> {
        let zero = Affine::<C>::zero();
        if pos >= B {
            return Err(AVError::WrongPosition(format!(
                "AV: Position has to be less than {}",
                B
            )));
        }
        if msg == zero {
            return Err(AVError::WrongMsg(String::from("AV: Message can't be zero")));
        }
        // assert that position is not already used
        match self.state {
            OracleState::Round1Ongoing => {
                if self.first_msgs[pos] != zero {
                    return Err(AVError::MsgAlreadySent(pos));
                }
            }
            OracleState::Round2Ongoing => {
                if self.second_msgs[pos] != zero {
                    return Err(AVError::MsgAlreadySent(pos));
                }
            }
            OracleState::Round1Completed => {
                return Err(AVError::WrongState(String::from(
                    "AV: Message can't be sent in Round1Completed state",
                )))
            }
            OracleState::Round2Completed => {
                return Err(AVError::WrongState(String::from(
                    "AV: Message can't be sent in Round2Completed state",
                )))
            }
            OracleState::Completed => {
                return Err(AVError::WrongState(String::from(
                    "AV: Message can't be sent in Completed state",
                )))
            }
        }

        Ok(())
    }

    pub fn register_msg(&mut self, msg: Affine::<C>, pos: usize) -> Result<(), AVError> {
        self.msg_validity(msg, pos)?;
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
////////////
/////////////////////////////////////////////////////////////////////////////
    pub fn output_first_round(&mut self) -> Vec<Affine::<C>> {
        assert_eq!(self.state, OracleState::Round1Completed);

        let mut x: Vec<C::ScalarField> = (0..B)
            .map(|i| if i == B - 1 { C::ScalarField::zero() }
                     else { C::ScalarField::one() })
        .collect();

        let mut outputs = vec![Affine::<C>::zero(); B];

        let cfg = MSMConfig::default();
        let mut projective_output = vec![Projective::<C>::zero(); 1];
        
        msm::msm(
            HostSlice::from_slice(&x),
            HostSlice::from_slice(&self.first_msgs),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();

        let mut affine_output = Affine::<C>::zero();
        C::to_affine(&projective_output[0], &mut affine_output);

        outputs[B-1] = affine_output;

        /*
           0 -1 -1 -1
           1  0 -1 -1
           1  1  0 -1
           1  1  1  0
        */
        for i in 0..(B - 1) {
            let idx = B - 2 - i;
            let projective_output = outputs[idx + 1].to_projective() - self.first_msgs[idx + 1].to_projective() - self.first_msgs[idx].to_projective();
            let mut affine_output = Affine::<C>::zero();
            C::to_affine(&projective_output, &mut affine_output);
            outputs[idx] = affine_output;
        }

        self.state = OracleState::Round2Ongoing;

        outputs
    }

    pub fn output_second_round(&mut self) -> Affine::<C> {
        assert_eq!(self.state, OracleState::Round2Completed);
        self.state = OracleState::Completed;

        let ones = vec![C::ScalarField::one(); B];

        let cfg = MSMConfig::default();
        let mut projective_output = vec![Projective::<C>::zero(); 1];
        
        msm::msm(
            HostSlice::from_slice(&ones),
            HostSlice::from_slice(&self.second_msgs),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();

        let mut affine_output = Affine::<C>::zero();
        C::to_affine(&projective_output[0], &mut affine_output);
        affine_output
    }
}

#[cfg(test)]
mod av_oracle_tests {
    use std::ops::Mul;

    use super::AVOracle;
    const B: usize = 1024;

    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg};
    use icicle_core::curve::{Curve,Affine};
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

        let round_outputs = av_oracle.output_first_round();
        let second_msgs: Vec<Affine::<Bn254CurveCfg>> = party_secrets
            .iter()
            .zip(round_outputs.iter())
            .map(|(mi, &ri)| Bn254CurveCfg::mul_scalar(ri.into(), *mi).into())
            .collect();

        for i in 0..B {
            av_oracle.register_msg(second_msgs[i], i).unwrap();
        }

        let output = av_oracle.output_second_round();
        assert_eq!(output, Affine::<Bn254CurveCfg>::zero());
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

        let round_outputs = av_oracle.output_first_round();

        // at least one party picks a new secret
        party_secrets[0] = ScalarCfg::generate_random(1)[0];

        let second_msgs: Vec<Affine::<Bn254CurveCfg>> = party_secrets
            .iter()
            .zip(round_outputs.iter())
            .map(|(mi, &ri)| Bn254CurveCfg::mul_scalar(ri.into(), *mi).into())
            .collect();

        for i in 0..B {
            av_oracle.register_msg(second_msgs[i], i).unwrap();
        }

        let output = av_oracle.output_second_round();
        assert_ne!(output, Affine::<Bn254CurveCfg>::zero());
    }
}