use icicle_core::curve::{Curve};
use icicle_core::pairing::Pairing;
use icicle_core::traits::{FieldImpl, Arithmetic};
use std::marker::PhantomData;

use crate::kzg::{Kzg, PK as KzgPk, VK as KzgVk};

use self::{
    structs::{Error, Instance, Proof, Witness},
    tr::Transcript,
};

use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{NTTDomain, NTTInitDomainConfig};
use std::ops::Mul;

use icicle_runtime::memory::HostSlice;
use icicle_core::ntt;

use ark_bn254::{Bn254, Fr as F, G1Projective, G2Projective};
use ark_ec::pairing::Pairing as OtherPairing;
use ark_ff::{FftField, Field, One};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain,
    Polynomial,
};

pub mod structs;
mod tr;

pub struct Argument<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    _e: PhantomData<(C1, C2, F)>,
}

impl<C1, C2, F> Argument<C1, C2, F>
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub fn prove<P>(
        instance: &Instance<C1>,
        witness: &Witness<C1, P>,
        pk: &KzgPk<C1, C2, F>,
    ) -> Proof<C1> 
    where
        P: UnivariatePolynomial<Field = C1::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic + Mul<Output = <C1 as Curve>::ScalarField>,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        let mut tr = Transcript::new(b"acc-transcript");
        
        let n: usize = instance.n;
        let omega = ntt::get_root_of_unity::<C1::ScalarField>(n.try_into().unwrap());

        tr.send_instance(instance);

        let p = C1::ScalarField::from_u32(5u32);
        debug_assert!((0..n).all(|i| p != omega.pow(i)), "p must not be in the domain");

        let acc_sh = witness.acc.eval(&(p * omega));
        let acc = witness.acc.eval(&p);
        let zh_eval = p.pow(n) - C1::ScalarField::one();
        
        let q = {
            let lhs = (acc_sh - instance.mu * acc) * (p - omega.pow(n - 1));
            lhs * zh_eval.inv()
        };
        
        tr.send_q(&q);
        let beta = tr.get_beta();

        let coeffs_0 = [C1::ScalarField::zero() - C1::ScalarField::one(), C1::ScalarField::one()];
        let divisor_poly_0 = P::from_coeffs(HostSlice::from_slice(&coeffs_0), 2);
        let (q_0, _) = witness.acc.divide(&divisor_poly_0);
        let q_0 = Kzg::commit(&pk, &q_0);
        
        let coeffs_1 = [C1::ScalarField::zero() - beta, C1::ScalarField::one()];
        let divisor_poly_1 = P::from_coeffs(HostSlice::from_slice(&coeffs_1), 2);
        let (q_1, r_1) = witness.acc.divide(&divisor_poly_1);
        let q_1 = Kzg::commit(&pk, &q_1);

        let coeffs_2 = [C1::ScalarField::zero() - (beta * omega), C1::ScalarField::one()];
        let divisor_poly_2 = P::from_coeffs(HostSlice::from_slice(&coeffs_2), 2);
        let (q_2, _) = witness.acc.divide(&divisor_poly_2);
        let q_2 = Kzg::commit(&pk, &q_2);
        
        Proof {
            q,
            acc_opening: witness.acc.eval(&beta),
            acc_shifted_opening: witness.acc.eval(&(beta * omega)),
            q_0,
            q_1,
            q_2,
        }
    }

    pub fn verify(
        instance: &Instance<C1>,
        proof: &Proof<C1>,
        vk: &KzgVk<C1, C2, F>,
    ) -> Result<(), Error> 
    where 
        <C1 as Curve>::ScalarField: Arithmetic + Mul<Output = <C1 as Curve>::ScalarField>,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        let n: usize = instance.n;
        let omega = ntt::get_root_of_unity::<C1::ScalarField>(n.try_into().unwrap());
        let mut tr = Transcript::new(b"acc-transcript");
        tr.send_instance(instance);
        tr.send_q(&proof.q);
        let beta = tr.get_beta();
        
        let kzg_at_one = Kzg::verify(
            &[instance.acc_cm],
            &[C1::ScalarField::one()],
            proof.q_0,
            C1::ScalarField::one(),
            C1::ScalarField::one(),
            vk,
        );
        
        if !kzg_at_one.is_ok() {
            return Err(Error::OpeningFailed);
        }

        let kzg_at_beta = Kzg::verify(
            &[instance.acc_cm],
            &[proof.acc_opening],
            proof.q_1,
            beta,
            C1::ScalarField::one(),
            vk,
        );

        if !kzg_at_beta.is_ok() {
            return Err(Error::OpeningFailed);
        }

        let kzg_at_beta_sh = Kzg::verify(
            &[instance.acc_cm],
            &[proof.acc_shifted_opening],
            proof.q_2,
            beta * omega,
            C1::ScalarField::one(),
            vk,
        );

        if !kzg_at_beta_sh.is_ok() {
            return Err(Error::OpeningFailed);
        }

        let zh_eval = beta.pow(n) - C1::ScalarField::one();
        let eq = {
            (proof.acc_shifted_opening - instance.mu * proof.acc_opening)
                * (beta - omega.pow(instance.n - 1))
                == proof.q * zh_eval
        };

        if !eq {
            return Err(Error::RelationCheck);
        }

        Ok(())
    }
}

#[cfg(test)]
mod acc_tests {
    use crate::{
        kzg::{Kzg, PK, VK},
        utils::srs::unsafe_setup_from_tau,
    };

    use super::{
        structs::{Instance, Witness},
        Argument,
    };
    
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    
    use icicle_runtime::memory::HostSlice;

    use icicle_core::{ntt::{
        self,
        NTTConfig, NTTDir, NTTInitDomainConfig,
    }};

    use icicle_core::ntt::ntt_inplace;
    use icicle_core::traits::FieldImpl;
    use icicle_core::curve::Curve;
    use icicle_core::ntt::NTTDomain;
    use icicle_bn254::curve::ScalarCfg;

    use std::marker::PhantomData;

    use icicle_bn254::polynomials::DensePolynomial as Bn254DensePolynomial;
    use icicle_core::polynomials::UnivariatePolynomial;

    use icicle_core::traits::Arithmetic;
    #[test]
    fn test_acc() {
        let n: usize = 16;
        let rou = ntt::get_root_of_unity::<Bn254ScalarField>(n.try_into().unwrap());

        ntt::initialize_domain::<Bn254ScalarField>(
            rou,
            &NTTInitDomainConfig::default(),
        )
        .unwrap();

        let tau = Bn254ScalarField::from_u32(17u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(n - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), _e: PhantomData, };
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());
        
        let mu = Bn254ScalarField::from_u32(100);
        let mut evals: Vec<Bn254ScalarField> = (0..n)
            .scan(Bn254ScalarField::one(), |state, _| {
                let out = *state;
                *state = *state * mu;
                Some(out)
            })
            .collect();

        let mut cfg = NTTConfig::<Bn254ScalarField>::default();
        ntt_inplace(HostSlice::from_mut_slice(&mut evals), NTTDir::kInverse, &cfg).unwrap();

        let acc_poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&evals), n);
        let acc_cm   = Kzg::commit(&pk, &acc_poly);

        let instance = Instance::<Bn254CurveCfg> { n, mu, acc_cm };
        let witness  = Witness { acc: acc_poly, _e: PhantomData,};

        let proof  = Argument::prove(&instance, &witness, &pk);
        let check  = Argument::verify(&instance, &proof, &vk);
        assert!(check.is_ok());
    }
}