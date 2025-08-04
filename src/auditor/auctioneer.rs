use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;

#[derive(Debug)]
pub enum Error {
    FirstRoundAuditFail,
    SecondRoundAuditFail,
    AuditFail { bid_at: usize, error_detail: Box<Error> },
}

pub struct AVAuditor<const B: usize, C: Curve + icicle_core::msm::MSM<C>> {
    pub(crate) first_round_msgs: [Affine::<C>; B],
    pub(crate) first_round_output: [Affine::<C>; B],
    pub(crate) second_round_msgs: [Affine::<C>; B],
    pub(crate) second_round_output: Affine::<C>,
}

impl<const B: usize, C: Curve + icicle_core::msm::MSM<C>> AVAuditor<B, C> {
    pub fn audit(&self) -> Result<(), Error> 
    where
        C: Curve<ScalarField: Arithmetic>
    { 
        // x = [0, -1, -1, ..., -1]
        let x: Vec<C::ScalarField> = (0..B)
            .map(|i| if i == 0 { C::ScalarField::one() }
                     else { C::ScalarField::zero() - C::ScalarField::one() })
        .collect();

        let mut outputs = vec![Affine::<C>::zero(); B];

        let cfg = MSMConfig::default();
        let mut projective_output = vec![Projective::<C>::zero(); 1];
        
        msm::msm(
            HostSlice::from_slice(&x),
            HostSlice::from_slice(&self.first_round_msgs),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();

        let mut affine_output = Affine::<C>::zero();
        C::to_affine(&projective_output[0], &mut affine_output);

        outputs[0] = affine_output;

        /*
           0 -1 -1 -1
           1  0 -1 -1
           1  1  0 -1
           1  1  1  0
        */
        for i in 1..B {
            let projective_output = outputs[i - 1].to_projective() + self.first_round_msgs[i - 1].to_projective() + self.first_round_msgs[i].to_projective();
            let mut affine_output = Affine::<C>::zero();
            C::to_affine(&projective_output, &mut affine_output);
            outputs[i] = affine_output;
        }

        let first_round_result = outputs;

        if first_round_result != self.first_round_output {
            return Err(Error::FirstRoundAuditFail);
        }

        let ones = vec![C::ScalarField::one(); B];
        
        let cfg = MSMConfig::default();
        let mut projective_output = vec![Projective::<C>::zero(); 1];
        
        msm::msm(
            HostSlice::from_slice(&ones),
            HostSlice::from_slice(&self.second_round_msgs),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();

        let mut affine_output = Affine::<C>::zero();
        C::to_affine(&projective_output[0], &mut affine_output);

        let second_round_result = affine_output;

        if second_round_result != self.second_round_output {
            return Err(Error::SecondRoundAuditFail);
        }

        Ok(())
    }
}

pub struct AuctioneerAuditor<const N: usize, const B: usize, C: Curve  + icicle_core::msm::MSM<C>> {
    av_auditor: [AVAuditor<B, C>; N],
}

impl<const N: usize, const B: usize, C: Curve + icicle_core::msm::MSM<C>> AuctioneerAuditor<N, B, C> {
    pub fn audit(&self) -> Result<(), Error> 
    where
        C: Curve<ScalarField: Arithmetic>
    {
        for i in 0..N {
            match self.av_auditor[i].audit() {
                Ok(()) => continue,
                Err(inner_error) => {
                    return Err(Error::AuditFail {
                        bid_at: i,
                        error_detail: Box::new(inner_error),
                    });
                }
            }
        }

        Ok(())
    }
}