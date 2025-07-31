use icicle_core::curve::{Curve,Affine};
use icicle_core::traits::FieldImpl;
use icicle_core::pairing::Pairing;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_runtime::memory::HostSlice;
use crate::utils::get_coeffs_of_poly;
use icicle_core::traits::Arithmetic;
use icicle_core::ntt::{get_root_of_unity,initialize_domain,NTTInitDomainConfig,NTTDomain,release_domain};

use structs::{Instance, Proof, Witness};

use crate::kzg::{Kzg, PK as KzgPk, VK as KzgVk};
use crate::utils::is_pow_2;

use self::structs::Error;
use self::tr::Transcript;

pub mod structs;
mod tr;

pub struct UnivariateSumcheck<C1, C2, F, U>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    _e: PhantomData<(C1, C2, F, U)>,
}

impl<C1, C2, F, U> UnivariateSumcheck<C1, C2, F, U>
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    pub fn new_instance(
        n: usize,
        a_cm: Affine::<C1>,
        b_cm: Affine::<C1>,
        sum: C1::ScalarField,
    ) -> Instance<C1> {
        assert!(is_pow_2(n));
        Instance { n, a_cm, b_cm, sum }
    }

    pub fn prove(
        witness: &Witness<C1, U>,
        instance: &Instance<C1>,
        pk: &KzgPk<C1, C2, F>,
    ) -> Proof<C1> 
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        assert!(is_pow_2(instance.n));
        let mut tr = Transcript::new(b"univariate-sumcheck");

        tr.send_instance(instance);

        //let len = instance.n + 1;
        //let poly_domain = get_root_of_unity::<C1::ScalarField>(len as u64);
        //initialize_domain(poly_domain, &NTTInitDomainConfig::default()).unwrap();
        let ab = witness.a_poly.mul(&witness.b_poly);

        let len = instance.n + 1;
        let mut vanishing_poly_coeffs = Vec::with_capacity(len);

        vanishing_poly_coeffs.push(C1::ScalarField::zero() - C1::ScalarField::one());
        for _ in 0..(len - 2) {
            vanishing_poly_coeffs.push(C1::ScalarField::zero());
        }
        vanishing_poly_coeffs.push(C1::ScalarField::one());

        let domain_vanishing_poly = U::from_coeffs(HostSlice::from_slice(&vanishing_poly_coeffs), len);
        
        let (q, r) = ab.divide(&domain_vanishing_poly);
        
        let mut r_coeffs = get_coeffs_of_poly(&r);
        
        while r_coeffs.len() > 1 && r_coeffs.last() == Some(&C1::ScalarField::zero()) {
            r_coeffs.pop();
        }
        
        let r_mod_x = U::from_coeffs(HostSlice::from_slice(&r_coeffs[1..]), r_coeffs[1..].len());

        let r_cm = Kzg::commit(pk, &r_mod_x);
        let q_cm = Kzg::commit(pk, &q);
        tr.send_r_and_q(&r_cm, &q_cm);

        let opening_challenge = tr.get_opening_challenge();

        let a_opening = witness.a_poly.eval(&opening_challenge);
        let b_opening = witness.b_poly.eval(&opening_challenge);

        let r_opening = r_mod_x.eval(&opening_challenge);
        let q_opening = q.eval(&opening_challenge);
        
        //release_domain::<C1::ScalarField>().unwrap();

        tr.send_openings(&a_opening, &b_opening, &r_opening, &q_opening);

        let separation_challenge = tr.get_separation_challenge();
        let pi = Kzg::open(
            pk,
            &[witness.a_poly.clone(), witness.b_poly.clone(), r_mod_x, q],
            opening_challenge,
            separation_challenge,
        );

        Proof {
            r_cm,
            q_cm,
            a_opening,
            b_opening,
            r_opening,
            q_opening,
            batch_opening_proof: pi,
        }
    }

    pub fn verify(
        proof: &Proof<C1>,
        instance: &Instance<C1>,
        vk: &KzgVk<C1, C2, F>,
    ) -> Result<(), Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        assert!(is_pow_2(instance.n));

        let commitments = [
            instance.a_cm.clone(),
            instance.b_cm.clone(),
            proof.r_cm.clone(),
            proof.q_cm.clone(),
        ];

        let evaluations = [
            proof.a_opening,
            proof.b_opening,
            proof.r_opening,
            proof.q_opening,
        ];

        let mut tr = Transcript::new(b"univariate-sumcheck");

        tr.send_instance(instance);
        tr.send_r_and_q(&proof.r_cm, &proof.q_cm);
        let opening_challenge = tr.get_opening_challenge();
        tr.send_openings(
            &proof.a_opening,
            &proof.b_opening,
            &proof.r_opening,
            &proof.q_opening,
        );
        let separation_challenge = tr.get_separation_challenge();

        // check a, b, r, q kzg opening proofs
        let opening_result = Kzg::verify(
            &commitments,
            &evaluations,
            proof.batch_opening_proof,
            opening_challenge,
            separation_challenge,
            vk,
        );

        assert!(opening_result.is_ok());

        let lhs = proof.a_opening * proof.b_opening;

        // check sumcheck relation
        let rhs = {
            let n_inv = C1::ScalarField::from_u32((instance.n).try_into().unwrap()).inv();
            opening_challenge * proof.r_opening
                + instance.sum * n_inv
                + proof.q_opening * (opening_challenge.pow(instance.n) - C1::ScalarField::one())
        };

        assert_eq!(lhs, rhs); //sometimes may fail

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;

    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_runtime::memory::HostSlice;
    
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_core::curve::Curve;

    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;
    
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity};

    use crate::utils::srs::unsafe_setup_from_tau;
    use std::marker::PhantomData;

    use super::*;
    use crate::kzg::{PK, VK};

    use icicle_core::ntt::{initialize_domain, release_domain};

    #[test]
    fn sumcheck_test() {
        let n = 16;

        let domain = get_root_of_unity::<Bn254ScalarField>((n * n).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

        //let domain = get_root_of_unity::<Bn254ScalarField>(n);

        let a_coeffs = ScalarCfg::generate_random(n as usize);
        let a_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&a_coeffs), n as usize);

        let cfg = NTTConfig::<Bn254ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut a_evals = vec![Bn254ScalarField::zero(); n as usize];

        ntt(
            HostSlice::from_slice(&a_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut a_evals),
        )
        .unwrap();

        let b_coeffs = ScalarCfg::generate_random(n as usize);
        let b_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&b_coeffs), n as usize);

        let mut b_evals = vec![Bn254ScalarField::zero(); n as usize];
        //domain
        ntt(
            HostSlice::from_slice(&b_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut b_evals),
        )
        .unwrap();

        //release_domain::<Bn254ScalarField>().unwrap();

        let mut sum = Bn254ScalarField::zero();
        for i in 0..a_evals.len() {
            sum = sum + (a_evals[i] * b_evals[i]);
        }

        let witness = Witness {
            a_poly: a_poly.clone(),
            b_poly: b_poly.clone(),
            e: PhantomData,
        };

        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(n as usize - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };

        let a_cm = Kzg::commit(&pk, &a_poly);
        let b_cm = Kzg::commit(&pk, &b_poly);

        let instance = Instance::<Bn254CurveCfg> { n: n as usize, a_cm, b_cm, sum };
        
        let proof = UnivariateSumcheck::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::prove(&witness, &instance, &pk);
        
        let result = UnivariateSumcheck::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&proof, &instance, &vk);
        
        release_domain::<Bn254ScalarField>().unwrap();

        assert!(result.is_ok());
    }
}