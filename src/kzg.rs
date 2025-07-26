use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;

use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;

use icicle_core::polynomials::UnivariatePolynomial;
use crate::utils::get_coeffs_of_poly;

use icicle_core::traits::Arithmetic;

use std::collections::BTreeMap;

#[derive(Debug)]
pub enum Error {
    PairingNot0,
}

pub struct Kzg<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub _e: PhantomData<(C1, C2, F)>,
}

pub struct PK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub srs: Vec<Affine<C1>>,
    pub _e: PhantomData<(C2, F)>,
}

pub struct VK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub g2: Affine<C2>,
    pub x_g2: Affine<C2>,
    _e: PhantomData<(C1, F)>,
}

pub struct DegreeCheckVK<C1, C2, F> 
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub pk_max_degree: usize,
    pub shifts: BTreeMap<usize, Affine::<C2>>,
    pub _e: PhantomData<(C1, F)>,
}

impl<C1, C2, F> DegreeCheckVK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub fn get_shift(&self, degree_bound: usize) -> Option<&Affine::<C2>> {
        let shift_factor = self.pk_max_degree - degree_bound;
        self.shifts.get(&shift_factor)
    }
}

impl<C1, C2, F> VK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub fn new(x_g2: Affine<C2>) -> Self {
        Self {
            g2: C2::get_generator().into(),
            x_g2,
            _e: PhantomData,
        }
    }
}

impl<C1, C2, F> Kzg<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    C1: icicle_core::msm::MSM<C1>,
{
    pub fn commit<P>(
        pk: &PK<C1, C2, F>,
        poly: &P
    ) -> Affine<C1> 
    where
        P: UnivariatePolynomial<Field = C1::ScalarField>,
    {
        if pk.srs.len() - 1 < poly.degree().try_into().unwrap() {
            panic!(
                "SRS size too small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                pk.srs.len()
            );
        }

        let mut cfg = MSMConfig::default();
        let mut projective_output = vec![Projective::<C1>::zero(); 1];
        let coeffs = get_coeffs_of_poly(poly);

        msm::msm(
            HostSlice::from_slice(&coeffs),
            HostSlice::from_slice(&pk.srs[..coeffs.len()]),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();
        
        let mut affine_output = Affine::<C1>::zero();
        C1::to_affine(&projective_output[0], &mut affine_output);
        affine_output
    }

    pub fn open<P>(
        pk: &PK<C1, C2, F>,
        polys: &[P],
        opening_challenge: C1::ScalarField,
        separation_challenge: C1::ScalarField,
    ) -> Affine<C1>
    where
        P: UnivariatePolynomial<Field = C1::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        let powers_of_gamma = std::iter::successors(Some(separation_challenge), |p| {
            Some(*p * separation_challenge)
        });

        let mut batched = polys[0].clone();
        for (p_i, gamma_pow_i) in polys.iter().skip(1).zip(powers_of_gamma) {
            let gamma_poly = p_i.mul_by_scalar(&gamma_pow_i);
            batched = gamma_poly.add(&batched); //&gamma_poly + &batched;
        }
        
        let coeffs = [C1::ScalarField::zero() - opening_challenge, C1::ScalarField::one()];
        let divisor_poly = P::from_coeffs(HostSlice::from_slice(&coeffs), 2);

        let (q, _) = batched.divide(&divisor_poly);

        if pk.srs.len() - 1 < q.degree().try_into().unwrap() {
            panic!(
                "Batch open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                pk.srs.len()
            );
        }

        Kzg::commit(pk, &q)
    }

    pub fn verify(
        commitments: &[Affine<C1>],
        evaluations: &[C1::ScalarField],
        opening_proof: Affine<C1>,
        opening_challenge: C1::ScalarField,
        separation_challenge: C1::ScalarField,
        vk: &VK<C1, C2, F>,
    ) -> Result<(), Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        assert_eq!(commitments.len(), evaluations.len());
        let powers_of_gamma: Vec<_> = std::iter::successors(Some(C1::ScalarField::one()), |p| {
            Some(*p * separation_challenge)
        })
        .take(commitments.len())
        .collect();
        
        let mut cfg = MSMConfig::default();
        cfg.is_async = false;
        let mut batched_commitment = vec![Projective::<C1>::zero(); 1];
        msm::msm(
            HostSlice::from_slice(&powers_of_gamma),
            HostSlice::from_slice(commitments),
            &cfg,
            HostSlice::from_mut_slice(&mut batched_commitment),
        )
        .unwrap();
        
        let batched_eval: C1::ScalarField = evaluations
            .iter()
            .zip(powers_of_gamma.iter())
            .map(|(&ei, &gamma_i)| ei * gamma_i)
            .fold(C1::ScalarField::zero(), |acc, x| acc + x);
        
        /*
            (p(X) - y) = q(X)(X - z)
            p(X) - y = q(X)•X - q(X)z
            p(X) - y + q(X)z = q(X)•X
            e([p] + z[q] - y[1], [1]) = e([q], [x])
            e([p] + z[q] - y[1], [1])*e([q], -[x]) = 0
        */
        
        let opening = opening_proof.to_projective() * opening_challenge;
        let minus_eval = C1::ScalarField::zero() - batched_eval;
        let y_root = C1::get_generator() * minus_eval;

        let lhs = batched_commitment[0] + opening + y_root;
        let lhs_affine: Affine<C1> = lhs.into();
        
        let e1: F = C1::pairing(&lhs_affine, &vk.g2).expect("pairing failed");
        let e2: F = C1::pairing(&opening_proof, &vk.x_g2).expect("pairing failed");

        if e1 != e2 {
            return Err(Error::PairingNot0);
        }

        Ok(())
    }
}

#[cfg(test)]
mod test_kzg {
    use icicle_bn254::polynomials::DensePolynomial as Bn254DensePolynomial;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;

    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_runtime::memory::HostSlice;
    
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg, G1Projective, G1Affine, G2Projective};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_core::curve::{Curve,Affine,Projective};

    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;
    
    use icicle_core::{msm, msm::MSMConfig};
    
    use crate::utils::srs::unsafe_setup_from_tau;
    use crate::utils::get_coeffs_of_poly;

    use super::{Kzg, PK, VK};
    use std::marker::PhantomData;

    use icicle_core::traits::Arithmetic;
    use std::collections::BTreeMap;
    use icicle_core::pairing::Pairing;
    
    #[test]
    fn test_kzg_real() {
        let size = 16;
        let coeffs = ScalarCfg::generate_random(size);
        let poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&coeffs), size);
        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(size - 1, tau);
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), _e: PhantomData, };

        let mut cfg = MSMConfig::default();
        let mut output = vec![G1Projective::zero(); size];
        
        let commit = Kzg::commit(&pk, &poly);

        let z = Bn254ScalarField::from_u32(10u32);
        let gamma = Bn254ScalarField::from_u32(20u32);
        
        let polys = [poly.clone()];
        let q = Kzg::open(&pk, &polys, z, gamma);

        let eval = poly.eval(&z);
        
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());
        
        let verify_result = Kzg::verify(&[commit], &[eval], q, z, gamma, &vk);
        assert!(verify_result.is_ok());
    }
    
    #[test]
    fn test_batched_kzg_real() {
        let size = 16;
        let a_coeffs = ScalarCfg::generate_random(size);
        let a_poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&a_coeffs), size);

        let b_coeffs = ScalarCfg::generate_random(size);
        let b_poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&b_coeffs), size);

        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(size - 1, tau);

        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), _e: PhantomData, };

        let mut cfg = MSMConfig::default();
        let mut output = vec![G1Projective::zero(); size];
        
        let a_commit = Kzg::commit(&pk, &a_poly);
        let b_commit = Kzg::commit(&pk, &b_poly);

        let z = Bn254ScalarField::from_u32(10u32);
        let gamma = Bn254ScalarField::from_u32(20u32);
        
        let polys = [a_poly.clone(), b_poly.clone()];
        let q = Kzg::open(&pk, &polys, z, gamma);

        let a_eval = a_poly.eval(&z);
        let b_eval = b_poly.eval(&z);

        let verify_result = Kzg::verify(&[a_commit, b_commit], &[a_eval, b_eval], q, z, gamma, &vk);
        assert!(verify_result.is_ok());
    }

    #[test]
    fn test_kzg_related_functions() {
        let size = 16;
        let coeffs = ScalarCfg::generate_random(size);
        let poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&coeffs), size);
        let tau = Bn254ScalarField::from_u32(100u32);

        let mut scalars = Vec::with_capacity(size);
        let mut acc = Bn254ScalarField::from_u32(1u32);
        for _ in 0..size {
            scalars.push(acc);
            acc = acc * tau;
        }
        assert_eq!(scalars.len(),size);
        
        let projective_g1: G1Projective = Bn254CurveCfg::get_generator();
        
        let mut srs = Vec::with_capacity(size);
        for i in &scalars {
            let projective_base = Bn254CurveCfg::mul_scalar(projective_g1, *i);
            let mut affine_base: G1Affine = G1Affine::zero();
            Bn254CurveCfg::to_affine(&projective_base, &mut affine_base);
            srs.push(affine_base);
        }
        
        let mut cfg = MSMConfig::default();
        let mut output = vec![G1Projective::zero(); size];
        
        msm::msm(
            HostSlice::from_slice(&coeffs),
            HostSlice::from_slice(&srs),
            &cfg,
            HostSlice::from_mut_slice(&mut output),
        )
        .unwrap();

        let w = Bn254ScalarField::from_u32(10u32);

        let eval = poly.eval(&w);
    }

    #[test]
    fn test_degree_bound() {
        let n = 16;
        let a_coeffs = ScalarCfg::generate_random(n);
        let a_poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&a_coeffs), n);

        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(2 * n - 1, tau);

        let shift_factor = srs.len() - 1 - (n - 1);

        let tau_pow_shift = Bn254G2CurveCfg::mul_scalar(Bn254G2CurveCfg::get_generator(), (tau.pow(shift_factor)));
        let mut degree_check_vk_map: BTreeMap<usize, Projective<Bn254G2CurveCfg>> = BTreeMap::new();
        degree_check_vk_map.insert(shift_factor, tau_pow_shift);

        // we want to check that a is of degree <= n-1
        let a_degree = {
            let mut coeffs = get_coeffs_of_poly(&a_poly);
            let mut shifted_coeffs = vec![Bn254ScalarField::zero(); shift_factor];
            shifted_coeffs.append(&mut coeffs);
            Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&shifted_coeffs), shifted_coeffs.len())
        };

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), _e: PhantomData, };
        let a_cm = Kzg::commit(&pk, &a_poly);
        let a_degree_cm = Kzg::commit(&pk, &a_degree);
        
        let mut affine_output = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(degree_check_vk_map.get(&shift_factor).unwrap(), &mut affine_output);

        let lhs = Bn254CurveCfg::pairing(&a_cm, &affine_output).expect("pairing failed");

        let mut affine_output = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&Bn254G2CurveCfg::get_generator(), &mut affine_output);

        let rhs = Bn254CurveCfg::pairing(&a_degree_cm, &affine_output).expect("pairing failed");

        assert_eq!(lhs, rhs);
    }
}