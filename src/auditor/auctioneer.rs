use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;

pub struct AVAuditor<const B: usize, C: Curve + icicle_core::msm::MSM<C>> {
    pub(crate) first_round_msgs: [Affine::<C>; B],
    pub(crate) first_round_output: [Affine::<C>; B],
    pub(crate) second_round_msgs: [Affine::<C>; B],
    pub(crate) second_round_output: Affine::<C>,
}

impl<const B: usize, C: Curve + icicle_core::msm::MSM<C>> AVAuditor<B, C> {
    pub fn audit(&self) -> bool 
    where
        C: Curve<ScalarField: Arithmetic>
    {
        let mut x = [C::ScalarField::zero() - C::ScalarField::one(); B];
        x[0] = C::ScalarField::zero();

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

        let first_round_res = outputs;

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

        let second_round_res = affine_output;

        return first_round_res == self.first_round_output
            && second_round_res == self.second_round_output;
    }
}

pub struct AuctioneerAuditor<const N: usize, const B: usize, C: Curve  + icicle_core::msm::MSM<C>> {
    av_auditor: [AVAuditor<B, C>; N],
}

impl<const N: usize, const B: usize, C: Curve + icicle_core::msm::MSM<C>> AuctioneerAuditor<N, B, C> {
    pub fn audit(&self) -> bool where C: Curve<ScalarField: Arithmetic> {
        for i in 0..N {
            if !self.av_auditor[i].audit() {
                return false;
            }
        }

        return true;
    }
}