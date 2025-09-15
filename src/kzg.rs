//! KZG polynomial commitments over generic curves, with optional GPU accelerated MSMs.
//!
//! This module provides a KZG interface with:
//! - `commit` — commit to a univariate polynomial
//! - `open` — produce a single or batched opening proof at a point
//! - `verify` — verify one or more evaluations against commitments and a single opening proof

use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;

use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;
use icicle_core::polynomials::UnivariatePolynomial;
use crate::utils::{msm_gpu, my_msm, get_coeffs_of_poly, get_device_is_cpu_or_gpu};

use icicle_core::traits::Arithmetic;

use std::collections::BTreeMap;

use icicle_runtime::Device;

pub mod lagrange;

/// Errors that can occur while using the KZG API.
/// - `PairingCheckFailed` — the final pairing equation did not hold during verification.
/// - `SrsTooSmall { required, actual }` — the proving key's SRS does not support the requested polynomial degree.
#[derive(Debug)]
pub enum Error {
    /// Pairing equality failed in verification.
    PairingCheckFailed,
    /// The proving key’s SRS does not support the required polynomial degree.
    /// `required` is the maximum supported degree, `actual` is the requested degree.
    SrsTooSmall { required: usize, actual: usize },
}

/// Zero-sized type implementing the KZG API over a generic pairing-friendly curve.
///
/// The type parameters:
/// - `C1` — the G1 curve (commitments and proofs live here),
/// - `C2` — the G2 curve (verifier key lives here),
/// - `F` — the target field of the pairing.
pub struct Kzg<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub e: PhantomData<(C1, C2, F)>,
}

/// Proving key: structured reference string (SRS) in G1 used to commit and open.
/// Its length defines the supported maximum degree.
/// # Fields
///
/// - `srs`: Powers of the trapdoor in G1,
///
/// # Type parameters
///
/// - `C1`, `C2`, `F`: See [`Kzg`].
///
/// # Notes
///
/// The SRS length must be strictly larger than the maximum polynomial degree you want to commit.
///
pub struct PK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub srs: Vec<Affine<C1>>,
    pub e: PhantomData<(C2, F)>,
}

/// Verifying key for KZG: fixed elements in G2.
///
/// # Fields
///
/// - `g2`: The G2 generator `[1]_2`.
/// - `x_g2`: The trapdoor in G2, `[tau]_2`.
///
/// # Type parameters
///
/// - `C1`, `C2`, `F`: See [`Kzg`].
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

/// Verifier material for degree bound check in pi_pse.
///
/// Maintains precomputed **shift elements** in G2 so the verifier can check that a polynomial’s
/// degree is at most a bound.
///
/// # Fields
///
/// - `pk_max_degree`: Maximum degree supported by the proving key.
/// - `shifts`: Map from `shift_factor` to the corresponding G2 element.
///
/// # Type parameters
///
/// - `C1`, `C2`, `F`: See [`Kzg`].
pub struct DegreeCheckVK<C1, C2, F> 
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub pk_max_degree: usize,
    pub shifts: BTreeMap<usize, Affine::<C2>>,
    pub e: PhantomData<(C1, F)>,
}

impl<C1, C2, F> DegreeCheckVK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    /// Return the G2 **shift** for a given degree bound, if it exists.
    ///
    /// A polynomial with degree bound `d` uses shift factor `pk_max_degree - d`.
    ///
    /// # Parameters
    ///
    /// - `degree_bound`: The claimed upper bound on the polynomial’s degree.
    ///
    /// # Returns
    ///
    /// - `Some(&Affine<C2>)` if the corresponding shift is present.
    /// - `None` if no shift for that bound has been set up.
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
    /// Construct a verifying key from `[tau]_2`.
    ///
    /// # Parameters
    ///
    /// - `x_g2`: The element `[tau]_2`.
    ///
    /// # Returns
    ///
    /// A [`VK`] with `g2 = [1]_2` and `x_g2 = [tau]_2`.
    ///
    /// # Type parameters
    ///
    /// - `C1`, `C2`, `F`: See [`Kzg`].
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
    /// Commit to a univariate polynomial.
    ///
    /// Computes the KZG commitment ∑ a_i * [tau^i]_1, where `a_i` are coefficients of `poly`.
    /// Uses a multi-scalar multiplication (MSM) that may be offloaded to the GPU when available.
    ///
    /// # Parameters
    ///
    /// - `pk`: Proving key that holds the G1 SRS.
    /// - `poly`: Polynomial(s) to commit to.
    ///
    /// # Returns
    ///
    /// - `Ok(Affine<C1>)`: The commitment in G1.
    /// - `Err(Error::SrsTooSmall)`: If the SRS does not cover `deg(poly)`.
    ///
    /// # Type parameters
    ///
    /// - `P`: Polynomial type implementing [`UnivariatePolynomial`].
    pub fn commit<P>(
        pk: &PK<C1, C2, F>,
        poly: &P,
    ) -> Result<Affine<C1>, Error> 
    where
        P: UnivariatePolynomial<Field = C1::ScalarField>,
    {   
        let cpu_or_gpu = get_device_is_cpu_or_gpu();

        let poly_degree = poly.degree().try_into().expect("Polynomial degree conversion failed");
        let srs_len = pk.srs.len();

        if poly_degree > srs_len - 1 {
            return Err(Error::SrsTooSmall { required: srs_len - 1, actual: poly_degree });
        }

        let coeffs = get_coeffs_of_poly(poly);

        let projective_output = my_msm(&coeffs, &pk.srs[..coeffs.len()], cpu_or_gpu);
        
        let mut affine_output = Affine::<C1>::zero();
        C1::to_affine(&projective_output, &mut affine_output);
        Ok(affine_output)
    }

    /// Produce a opening proof at a challenge point.
    ///
    /// Given polynomials `p_0, …, p_{m-1}` and two challenges, z and gamma.
    /// This function folds the polynomials into one then computes an opening proof on the result polynomial.
    ///
    /// # Parameters
    ///
    /// - `pk`: Proving key.
    /// - `polys`: Polynomials to be opened together.
    /// - `opening_challenge`: The evaluation point `z`.
    /// - `batched_challenge`: The batching challenge `gamma`.
    ///
    /// # Returns
    ///
    /// - `Ok(Affine<C1>)`: The opening proof in G1.
    /// - `Err(Error::SrsTooSmall)`: If the SRS is too short for the derived quotient polynomial.
    ///
    /// # Type parameters
    ///
    /// - `P`: Polynomial type implementing [`UnivariatePolynomial`].
    pub fn open<P>(
        pk: &PK<C1, C2, F>,
        polys: &[P],
        opening_challenge: C1::ScalarField,
        batched_challenge: C1::ScalarField,
    ) -> Result<Affine<C1>, Error>
    where
        P: UnivariatePolynomial<Field = C1::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        let len_of_polys = polys.len();
        let mut powers_of_gamma = Vec::with_capacity(len_of_polys);
        // get gammas, gamma, gamma^2, ..., to fold many polys into one for batched opening
        if len_of_polys > 1 {
            let mut acc = batched_challenge;
            for _ in 0..len_of_polys {
                powers_of_gamma.push(acc);
                acc = acc * batched_challenge;
            }
        }
        
        let mut batched_poly = polys[0].clone();
        
        // random linear combination of polynomials
        if len_of_polys > 1 {
            for i in 1..len_of_polys {
                let gamma_poly = polys[i].mul_by_scalar(&powers_of_gamma[i - 1]);
                batched_poly = gamma_poly.add(&batched_poly);
            }
        }
        
        let coeffs = [C1::ScalarField::zero() - opening_challenge, C1::ScalarField::one()];
        let divisor_poly = P::from_coeffs(HostSlice::from_slice(&coeffs), 2);

        let (q, _) = batched_poly.divide(&divisor_poly);

        let q_degree = q.degree().try_into().expect("Polynomial degree conversion failed");
        let srs_len = pk.srs.len();

        if q_degree > srs_len - 1 {
            return Err(Error::SrsTooSmall { required: srs_len - 1, actual: q_degree });
        }

        Kzg::commit(pk, &q)
    }

    /// Verify one or more evaluations against commitments using a single opening proof.
    ///
    /// Batches `(commitment_i, evaluation_i)` pairs with `gamma`, forms the pairing equation
    ///
    /// e([C]_1 + z[Q]_1 - y[1]_1, [1]_2) = e([Q]_1, [tau]_2)
    ///
    /// where [Q]_1 is the opening proof and y is the batched evaluation.
    /// It is accepted iff the equality holds.
    ///
    /// # Parameters
    ///
    /// - `commitments`: Commitments in G1, one per evaluation.
    /// - `evaluations`: Claimed evaluations of polynomials.
    /// - `opening_proof`: The opening proof of KZG.
    /// - `opening_challenge`: Evaluation point `z`.
    /// - `batched_challenge`: Batching challenge `gamma`.
    /// - `vk`: Verifying key containing `[1]_2` and `[tau]_2`.
    ///
    /// # Returns
    ///
    /// - `Ok(())`: The proof verified successfully.
    /// - `Err(Error::PairingCheckFailed)`: The pairing check failed.
    pub fn verify(
        commitments: &[Affine<C1>],
        evaluations: &[C1::ScalarField],
        opening_proof: Affine<C1>,
        opening_challenge: C1::ScalarField,
        batched_challenge: C1::ScalarField,
        vk: &VK<C1, C2, F>,
    ) -> Result<(), Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        assert_eq!(commitments.len(), evaluations.len());
        let cpu_or_gpu = get_device_is_cpu_or_gpu();

        let len_of_commitments = commitments.len();
        let mut powers_of_gamma = Vec::with_capacity(commitments.len());
        
        // get gammas for batched verifing
        let mut acc = C1::ScalarField::one();
        for _ in 0..commitments.len() {
            powers_of_gamma.push(acc);
            acc = acc * batched_challenge;
        }

        let batched_commitment = my_msm(&powers_of_gamma, commitments, cpu_or_gpu);
        
        // random linear combination of evaluations
        let mut batched_eval = C1::ScalarField::zero();
        for i in 0..evaluations.len() {
            batched_eval = batched_eval + evaluations[i] * powers_of_gamma[i];
        }
        
        /*
            (p(X) - y) = q(X)(X - z)
            p(X) - y = q(X)•X - q(X)z
            p(X) - y + q(X)z = q(X)•X
            e([p] - y[1] + z[q], [1]) = e([q], [x])
        */
        
        // z[Q]
        let opening = opening_proof.to_projective() * opening_challenge;
        // -y
        let minus_eval = C1::ScalarField::zero() - batched_eval;
        // -y[1]_1
        let y_g = C1::get_generator() * minus_eval;
        
        // [C] - y[1] + z[Q] 
        let lhs = batched_commitment + y_g + opening;
        let lhs_affine: Affine<C1> = lhs.into();
        
        // e([C] - y[1] + z[Q], [1])
        let el: F = C1::pairing(&lhs_affine, &vk.g2).expect("pairing failed");
        // e([Q], [x])
        let er: F = C1::pairing(&opening_proof, &vk.x_g2).expect("pairing failed");

        if el != er {
            return Err(Error::PairingCheckFailed);
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
    
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg, G1Projective, G1Affine};
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
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        
        let commit = Kzg::commit(&pk, &poly).unwrap();

        let z = Bn254ScalarField::from_u32(10u32);
        let gamma = Bn254ScalarField::from_u32(20u32);
        
        let polys = [poly.clone()];
        let q = Kzg::open(&pk, &polys, z, gamma).unwrap();

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
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        
        let a_commit = Kzg::commit(&pk, &a_poly).unwrap();
        let b_commit = Kzg::commit(&pk, &b_poly).unwrap();

        let z = Bn254ScalarField::from_u32(10u32);
        let gamma = Bn254ScalarField::from_u32(20u32);
        
        let polys = [a_poly.clone(), b_poly.clone()];
        let q = Kzg::open(&pk, &polys, z, gamma).unwrap();

        let a_eval = a_poly.eval(&z);
        let b_eval = b_poly.eval(&z);

        let verify_result = Kzg::verify(&[a_commit, b_commit], &[a_eval, b_eval], q, z, gamma, &vk);
        assert!(verify_result.is_ok());
    }
    
    // it just is used to test related functions
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
        
        let cfg = MSMConfig::default();
        let mut output = vec![G1Projective::zero(); size];
        
        msm::msm(
            HostSlice::from_slice(&coeffs),
            HostSlice::from_slice(&srs),
            &cfg,
            HostSlice::from_mut_slice(&mut output),
        )
        .unwrap();

        let w = Bn254ScalarField::from_u32(10u32);

        let _eval = poly.eval(&w);
    }

    #[test]
    fn test_degree_bound() {
        let n = 16;
        let coeffs = ScalarCfg::generate_random(n);
        let poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&coeffs), n);

        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(2 * n - 1, tau);

        let shift_factor = srs.len() - 1 - (n - 1);

        let tau_pow_shift = Bn254G2CurveCfg::mul_scalar(Bn254G2CurveCfg::get_generator(), tau.pow(shift_factor));
        let mut degree_check_vk_map: BTreeMap<usize, Projective<Bn254G2CurveCfg>> = BTreeMap::new();
        degree_check_vk_map.insert(shift_factor, tau_pow_shift);

        // the degree bound, deg <= n - 1
        let degree = {
            let mut coeffs = get_coeffs_of_poly(&poly);
            let mut shifted_coeffs = vec![Bn254ScalarField::zero(); shift_factor];
            shifted_coeffs.append(&mut coeffs);
            Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&shifted_coeffs), shifted_coeffs.len())
        };

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        let commit = Kzg::commit(&pk, &poly).unwrap();
        let degree_commit = Kzg::commit(&pk, &degree).unwrap();
        
        let mut affine_output = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(degree_check_vk_map.get(&shift_factor).unwrap(), &mut affine_output);

        let lhs = Bn254CurveCfg::pairing(&commit, &affine_output).expect("pairing failed");

        let mut affine_output = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&Bn254G2CurveCfg::get_generator(), &mut affine_output);

        let rhs = Bn254CurveCfg::pairing(&degree_commit, &affine_output).expect("pairing failed");

        assert_eq!(lhs, rhs);
    }
}