use icicle_core::curve::{Curve,Affine};
use icicle_core::traits::FieldImpl;
use icicle_core::pairing::Pairing;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_runtime::memory::HostSlice;
use crate::utils::get_coeffs_of_poly;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, NTT, release_domain};
use icicle_core::traits::GenerateRandom;
use icicle_core::traits::Arithmetic;
use std::ops::{Mul, Add};

pub mod structs;
mod tr;

use self::structs::{Error, Instance, Proof, Witness};
use crate::{
    kzg::{DegreeCheckVK, Kzg, PK as KzgPk, VK as KzgVk},
    utils::{folding::compute_folding_coeffs, is_pow_2, evaluate_all_lagrange_coefficients},
    verifiable_folding_sumcheck::tr::Transcript,
};

pub struct Argument<C1, C2, F, U>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    _e: PhantomData<(C1, C2, F, U)>,
}

impl<C1, C2, F, U> Argument<C1, C2, F, U> 
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
    <U as UnivariatePolynomial>::FieldConfig: GenerateRandom<<U as UnivariatePolynomial>::Field>,
{
    pub fn sample_blinder(
        sum: C1::ScalarField,
        degree: usize,
        n: u64,
    ) -> U 
    where
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        let mut coeffs = <<U as UnivariatePolynomial>::FieldConfig as GenerateRandom<<U as UnivariatePolynomial>::Field>>::generate_random(degree + 1);

        let n_inv = C1::ScalarField::from_u32((n).try_into().unwrap()).inv();
        coeffs[0] = sum * n_inv;
        
        let blinder = U::from_coeffs(HostSlice::from_slice(&coeffs), coeffs.len());

        blinder
    }

    pub fn prove(
        instance: &Instance<C1>,
        witness: &Witness<C1, U>,
        pk: &KzgPk<C1, C2, F>,
    ) -> Proof<C1> 
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        assert!(is_pow_2(instance.n));

        let domain = get_root_of_unity::<C1::ScalarField>(instance.n.try_into().unwrap());
        let mut tr = Transcript::new(b"verifiable-folding-sumcheck");

        tr.send_instance(instance);

        let (b_1, b_2) = (<<U as UnivariatePolynomial>::FieldConfig as GenerateRandom<<U as UnivariatePolynomial>::Field>>::generate_random(1)[0], <<U as UnivariatePolynomial>::FieldConfig as GenerateRandom<<U as UnivariatePolynomial>::Field>>::generate_random(1)[0]);
        let s = C1::mul_scalar(instance.p_base.to_projective(), b_1) + C1::mul_scalar(instance.h_base.to_projective(), b_2);
        
        let mut s_affine = Affine::<C1>::zero();
        C1::to_affine(&s, &mut s_affine);

        let blinder = Self::sample_blinder(b_1, 1, instance.n as u64);
        let blinder_cm = Kzg::commit(pk, &blinder);
        tr.send_blinders(&s_affine, &blinder_cm);

        let c = tr.get_c();

        let z_1 = c * witness.x + b_1;
        let z_2 = c * witness.r + b_2;

        let b_evals = compute_folding_coeffs::<C1>(&instance.challenges);

        let mut cfg = NTTConfig::<C1::ScalarField>::default();
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut b_coeffs = vec![C1::ScalarField::zero(); b_evals.len()];
        ntt(
            HostSlice::from_slice(&b_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut b_coeffs),
        )
        .unwrap();
        release_domain::<C1::ScalarField>();
        let b = U::from_coeffs(HostSlice::from_slice(&b_coeffs), b_coeffs.len());
        
        // B(X) + ca(X)b(X)
        let a_nof_coeffs = witness.a.get_nof_coeffs();
        let b_nof_coeffs:u64 = b_coeffs.len().try_into().unwrap();
        assert_eq!(a_nof_coeffs, b_nof_coeffs);
        let max_domain_size:u64= a_nof_coeffs + b_nof_coeffs;
        let poly_domain = get_root_of_unity::<C1::ScalarField>(max_domain_size);
        initialize_domain(poly_domain, &NTTInitDomainConfig::default()).unwrap();

        let mut lhs_product = witness.a.mul(&b);
        lhs_product = lhs_product.mul_by_scalar(&c);
        let lhs: U = blinder.add(&lhs_product);
        
        let len = instance.n + 1;
        let mut vanishing_poly_coeffs = Vec::with_capacity(len);
        release_domain::<C1::ScalarField>();

        vanishing_poly_coeffs.push(C1::ScalarField::zero() - C1::ScalarField::one());
        for i in 0..(len - 2) {
            vanishing_poly_coeffs.push(C1::ScalarField::zero());
        }
        vanishing_poly_coeffs.push(C1::ScalarField::one());

        let domain_vanishing_poly = U::from_coeffs(HostSlice::from_slice(&vanishing_poly_coeffs), len);
        
        let poly_domain = get_root_of_unity::<C1::ScalarField>(len as u64);
        initialize_domain(poly_domain, &NTTInitDomainConfig::default()).unwrap();

        let (q, r) = lhs.divide(&domain_vanishing_poly);

        let mut r_coeffs = get_coeffs_of_poly(&r);

        while r_coeffs.len() > 1 && r_coeffs.last() == Some(&C1::ScalarField::zero()) {
            r_coeffs.pop();
        }

        assert_eq!(z_1, C1::ScalarField::from_u32(instance.n as u32) * r_coeffs[0]);

        let r_mod_x = U::from_coeffs(HostSlice::from_slice(&r_coeffs[1..]), r_coeffs[1..].len());
        
        // deg(r_mod_x) <= n - 2
        let r_degree = {
            let shift_factor = pk.srs.len() - 1 - (instance.n - 2);
            let mut coeffs = get_coeffs_of_poly(&r_mod_x);
            let mut shifted_coeffs = vec![C1::ScalarField::zero(); shift_factor];
            shifted_coeffs.append(&mut coeffs);
            while shifted_coeffs.len() > 1 && shifted_coeffs.last() == Some(&C1::ScalarField::zero()) {
                shifted_coeffs.pop();
            }
            U::from_coeffs(HostSlice::from_slice(&shifted_coeffs), shifted_coeffs.len())
        };

        let r_mod_x_cm = Kzg::commit(pk, &r_mod_x);
        let r_degree_cm = Kzg::commit(pk, &r_degree);
        let q_cm = Kzg::commit(pk, &q);

        tr.second_round(&z_1, &z_2, &r_mod_x_cm, &r_degree_cm, &q_cm);
        let opening_challenge = tr.get_opening_challenge();

        let a_opening = witness.a.eval(&opening_challenge);
        let blinder_opening = blinder.eval(&opening_challenge);

        let r_opening = r_mod_x.eval(&opening_challenge);
        let q_opening = q.eval(&opening_challenge);
        release_domain::<C1::ScalarField>();
        
        tr.send_openings(&a_opening, &blinder_opening, &r_opening, &q_opening);

        let separation_challenge = tr.get_separation_challenge();
        let batch_opening_proof = Kzg::open(
            pk,
            &[witness.a.clone(), blinder.clone(), r_mod_x.clone(), q.clone()],
            opening_challenge,
            separation_challenge,
        );

        Proof {
            // round 1
            s: s_affine,
            blinder_cm,

            // round 2
            z_1,
            z_2,
            r_cm: r_mod_x_cm,
            r_degree_cm,
            q_cm,

            // round 3
            a_opening,
            blinder_opening,
            r_opening,
            q_opening,

            // round 4
            batch_opening_proof,
        }
    }

    pub fn verify(
        instance: &Instance<C1>,
        proof: &Proof<C1>,
        vk: &KzgVk<C1, C2, F>,
        degree_check_vk: &DegreeCheckVK<C1, C2, F>,
    ) -> Result<(), Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        let domain = get_root_of_unity::<C1::ScalarField>(instance.n.try_into().unwrap());
        let mut tr = Transcript::new(b"verifiable-folding-sumcheck");

        tr.send_instance(instance);
        tr.send_blinders(&proof.s, &proof.blinder_cm);
        let c = tr.get_c();

        tr.second_round(
            &proof.z_1,
            &proof.z_2,
            &proof.r_cm,
            &proof.r_degree_cm,
            &proof.q_cm,
        );
        let opening_challenge = tr.get_opening_challenge();

        tr.send_openings(
            &proof.a_opening,
            &proof.blinder_opening,
            &proof.r_opening,
            &proof.q_opening,
        );
        let separation_challenge = tr.get_separation_challenge();

        let commitments = [
            instance.a_cm.clone(),
            proof.blinder_cm.clone(),
            proof.r_cm.clone(),
            proof.q_cm.clone(),
        ];

        let evaluations = [
            proof.a_opening,
            proof.blinder_opening,
            proof.r_opening,
            proof.q_opening,
        ];

        let kzg_check = Kzg::verify(
            &commitments,
            &evaluations,
            proof.batch_opening_proof,
            opening_challenge,
            separation_challenge,
            vk,
        );

        if !kzg_check.is_ok() {
            return Err(Error::OpeningFailed);
        }

        let p_c = C1::mul_scalar(instance.pedersen.to_projective(), c);
        let leq = p_c + proof.s.to_projective();

        let p_z_1 = C1::mul_scalar(instance.p_base.to_projective(), proof.z_1);
        let h_z_2 = C1::mul_scalar(instance.h_base.to_projective(), proof.z_2);
        let req = p_z_1 + h_z_2;
        
        if leq != req {
            return Err(Error::PedersenOpeningFailed);
        }

        let b_evals = compute_folding_coeffs::<C1>(&instance.challenges);

        let lagrange_evals = evaluate_all_lagrange_coefficients::<C1>(domain, opening_challenge, instance.n);
        
        let mut b_opening = C1::ScalarField::zero();
        for i in 0..b_evals.len() {
            let bi = b_evals[i];
            let pi = lagrange_evals[i];
            b_opening = b_opening + (bi * pi);
        }

        let lhs = proof.blinder_opening + c * proof.a_opening * b_opening;

        let rhs = {
            let n_inv = C1::ScalarField::from_u32(instance.n as u32).inv();
            opening_challenge * proof.r_opening
                + proof.z_1 * n_inv
                + proof.q_opening * (opening_challenge.pow(instance.n) - C1::ScalarField::one())
        };
        
        if lhs != rhs {
            return Err(Error::RelationCheckFailed);
        }

        let shift = degree_check_vk.get_shift(instance.n - 2);
        let shift = match shift {
            Some(value) => Ok(*value),
            None => Err(Error::DegreeCheckShiftMissing),
        }?;

        let e1: F = C1::pairing(&proof.r_cm, &shift).expect("pairing failed");
        let e2: F = C1::pairing(&proof.r_degree_cm, &vk.g2).expect("pairing failed");
        
        if e1 != e2 {
            return Err(Error::DegreeCheckFailed);
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod verifiable_folding_sumcheck_tests {
    use std::{collections::BTreeMap, ops::Mul};

    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::traits::FieldImpl;
    use std::marker::PhantomData;
    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, ntt_inplace, NTT, release_domain};
    use icicle_runtime::memory::HostSlice;
    use icicle_core::traits::Arithmetic;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;

    use crate::{
        kzg::{DegreeCheckVK, Kzg, PK, VK},
        utils::{folding::compute_folding_coeffs, srs::unsafe_setup_from_tau},
    };

    use super::{
        structs::{Instance, Witness},
        Argument,
    };

    #[test]

    fn test_sumcheck() {
        let log_n = 4;
        let n = 1 << log_n;
        let domain = get_root_of_unity::<Bn254ScalarField>(n);

        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>((n - 1) as usize, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let mut x_g2_affine = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&x_g2, &mut x_g2_affine);

        // we will be checking for R <= n - 2
        let shift_factor = srs.len() - 1 - ((n - 2) as usize);
        let tau_pow_shift = Bn254G2CurveCfg::get_generator() * tau.pow(shift_factor);

        let mut tau_pow_shift_affine = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&tau_pow_shift, &mut tau_pow_shift_affine);

        let mut degree_check_vk_map: BTreeMap<usize, Affine<Bn254G2CurveCfg>> = BTreeMap::new();
        degree_check_vk_map.insert(shift_factor, tau_pow_shift_affine);

        let degree_check_vk = DegreeCheckVK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> {
            pk_max_degree: srs.len() - 1,
            shifts: degree_check_vk_map,
            _e: PhantomData,
        };

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), _e: PhantomData, };
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2_affine);

        let a_coeffs = ScalarCfg::generate_random(n as usize);
        let a_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&a_coeffs), a_coeffs.len());

        let mut cfg = NTTConfig::<Bn254ScalarField>::default();
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut a_evals = vec![Bn254ScalarField::zero(); n as usize];

        ntt(
            HostSlice::from_slice(&a_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut a_evals),
        )
        .unwrap();

        let challenges = ScalarCfg::generate_random(log_n);

        let g = Bn254CurveCfg::get_generator();
        let p = Bn254ScalarField::from_u32(200u32);
        let h = Bn254ScalarField::from_u32(300u32);
        
        let p_base = Bn254CurveCfg::mul_scalar(g, p);
        let h_base = Bn254CurveCfg::mul_scalar(g, h);

        let r = Bn254ScalarField::from_u32(10u32);

        let b_evals = compute_folding_coeffs::<Bn254CurveCfg>(&challenges);

        let mut x = Bn254ScalarField::zero();
        let len = a_evals.len();

        for i in 0..len {
            let ai = a_evals[i];
            let bi = b_evals[i];

             x = x + (ai * bi);
        }

        let pedersen = Bn254CurveCfg::mul_scalar(p_base, x) + Bn254CurveCfg::mul_scalar(h_base, r);

        let a_cm = Kzg::commit(&pk, &a_poly);

        let mut p_base_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&p_base, &mut p_base_affine);
        let mut h_base_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&h_base, &mut h_base_affine);
        let mut pedersen_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&pedersen, &mut pedersen_affine);

        let instance = Instance::<Bn254CurveCfg> {
            n: n as usize,
            p_base: p_base_affine,
            h_base: h_base_affine,
            a_cm,
            pedersen: pedersen_affine,
            challenges,
        };

        let witness = Witness { a: a_poly, r, x };

        let proof = Argument::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::prove(&instance, &witness, &pk);
        let res = Argument::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&instance, &proof, &vk, &degree_check_vk);

        assert!(res.is_ok());
    }
}