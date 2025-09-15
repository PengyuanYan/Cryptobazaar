//! Zero-knowledge log-derivative check over a masked subset of an evaluation domain.
//!
//! We have a polynomial `f(X)` defined over an `N`-sized domain. The last `B` points are
//! *blinders*. We want to prove, in zero knowledge, that for a random challenge `β`:
//!
//! ```text
//!     Σ_{i=0}^{N-B-1} 1/(β + f_i)  =  R(β)
//! ```
//!
//! where `R` is a public (or verifier-computed) rational function encoding the expected
//! log-derivative relation. A selector `s(X)` masks the active rows: `s = 1` on the first
//! `N-B` points and `0` on the last `B`.
//!
//! The prover commits to witness polynomials, transcript samples challenges,
//! prover constructs helper polynomials `b(X)` and `q(X)` so that a single batched KZG
//! opening enforces the relation at a random evaluation point.
use icicle_core::curve::Curve;
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, NTT, release_domain};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;
use crate::kzg::{Kzg, PK as KzgPk, VK as KzgVk};
use crate::utils::{my_ntt, get_device_is_cpu_or_gpu, get_coeffs_of_poly};
use std::marker::PhantomData;

use self::{
    structs::{Error, Instance, Proof, ProverIndex, VerifierIndex, Witness},
    tr::Transcript,
};

pub mod structs;
mod tr;

pub struct Argument<const N: usize, const B: usize, C1, C2, F, U> 
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    _e: PhantomData<(C1, C2, F, U)>,
}

impl<const N: usize, const B: usize, C1, C2, F, U> Argument<N, B, C1, C2, F, U> 
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    /// Computes the selector polynomial `s(X)` with `1` on the first `N−B` rows and `0` on the
    /// last `B`, commits to it with KZG, and returns the commitment.
    pub fn index_v(pk: &KzgPk<C1,C2,F>) -> VerifierIndex<C1> 
    where
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        let cpu_or_gpu = get_device_is_cpu_or_gpu();
        let mut zeros = vec![C1::ScalarField::zero(); B];
        let mut s_evals = vec![C1::ScalarField::one(); N - B];

        s_evals.append(&mut zeros);
        
        //domain
        let coeffs = my_ntt::<C1>(&s_evals, C1::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);

        let s = U::from_coeffs(HostSlice::from_slice(&coeffs), N);
        let s_cm = Kzg::commit(pk, &s).unwrap();

        VerifierIndex { s_cm: s_cm.into() }
    }

    /// Precomputes `s(X)` in coefficient form and its evaluations on the coset domain via NTT,
    /// so the prover can reuse them across proofs.
    pub fn index_p() -> ProverIndex<C1, U> 
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {  
        let cpu_or_gpu = get_device_is_cpu_or_gpu();
        let g = C1::ScalarField::from_u32(5u32);

        let mut zeros = vec![C1::ScalarField::zero(); B];
        let mut s_evals = vec![C1::ScalarField::one(); N - B];

        s_evals.append(&mut zeros);

        //domain
        let mut s_coeffs = my_ntt::<C1>(&s_evals, C1::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let mut s_coeffs_clone = s_coeffs.clone();
        
        //domain
        let mut s_coset_evals = my_ntt::<C1>(&s_coeffs, g, NTTDir::kForward, cpu_or_gpu);

        ProverIndex {
            s: U::from_coeffs(HostSlice::from_slice(&s_coeffs_clone), N),
            s_coset_evals: s_coset_evals,
        }
    }

    /// Create a proof of the masked log-derivative identity.
    ///
    /// 1. Commit to witness polynomials (e.g., `f`).
    /// 2. Sample challenges (`β`, `μ`, separation challenge).
    /// 3. Build helper polynomials `b(X)` and `q(X)` tying `Σ 1/(β+f_i)` to the claimed `R(β)`
    ///    while excluding the last `B` (via `s`).
    /// 4. Open all polynomials at a random point and batch-verify with KZG.
    pub fn prove(
        index: &ProverIndex<C1, U>,
        v_index: &VerifierIndex<C1>,
        instance: &Instance<C1>,
        witness: &Witness<C1, U>,
        pk: &KzgPk<C1, C2, F>,
    ) -> Proof<C1> 
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        let cpu_or_gpu = get_device_is_cpu_or_gpu();
        let g = C1::ScalarField::from_u32(5u32);

        let mut f_coeffs = get_coeffs_of_poly(&witness.f);

        let cfg = NTTConfig::<C1::ScalarField>::default();

        //domain
        let mut f_evals = my_ntt::<C1>(&f_coeffs, C1::ScalarField::one(), NTTDir::kForward, cpu_or_gpu);

        //domain
        let mut f_coset_evals = my_ntt::<C1>(&f_coeffs, g, NTTDir::kForward, cpu_or_gpu);

        let mut tr = Transcript::new_transcript(b"log-derivative");
        tr.send_v_index(v_index);
        tr.send_instance(instance);

        let mut blinders = f_evals[N - B..].to_vec();
        for i in 0..blinders.len() {
            blinders[i] = blinders[i].inv();
        }

        let mut sum = C1::ScalarField::zero();
        for i in 0..blinders.len() {
            sum = sum + blinders[i];
        }

        tr.send_blinders_sum(&sum);
        let gamma = tr.get_beta();

        let mut b_evals = Vec::with_capacity(N - B);
        for i in 0..(N - B) {
            b_evals.push((f_evals[i] + gamma).inv());
        }
        b_evals.append(&mut blinders);
        
        //domain
        let mut b_coeffs = my_ntt::<C1>(&b_evals, C1::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);

        let b = U::from_coeffs(HostSlice::from_slice(&b_coeffs), N);
        let b_cm = Kzg::commit(pk, &b).unwrap();

        //domain
        let mut b_coset_evals = my_ntt::<C1>(&b_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        
        let x = C1::ScalarField::from_u32(5u32);
        let zh_coset_inv = (x.pow(N.try_into().unwrap()) - C1::ScalarField::one()).inv();

        let mut q_coset_evals = Vec::with_capacity(N);
        for i in 0..N {
            let b_i = b_coset_evals[i];
            let f_i = f_coset_evals[i];
            let s_i = index.s_coset_evals[i];

            let result = (b_i * (s_i * gamma + f_i) - C1::ScalarField::one()) * zh_coset_inv;
            q_coset_evals.push(result);
        }
        
        //domain
        let mut q_coeffs = my_ntt::<C1>(&q_coset_evals, g, NTTDir::kInverse, cpu_or_gpu);

        let q = U::from_coeffs(HostSlice::from_slice(&q_coeffs), N);
        let q_cm = Kzg::commit(pk, &q).unwrap();

        tr.send_b_and_q(&b_cm, &q_cm);
        let mu = tr.get_mu();

        let f_opening = witness.f.eval(&mu);
        let s_opening = index.s.eval(&mu);
        let b_opening = b.eval(&mu);
        let q_opening = q.eval(&mu);

        tr.send_openings(&f_opening, &s_opening, &b_opening, &q_opening);
        let separation_challenge = tr.get_separation_challenge();

        let coeffs = [C1::ScalarField::zero(), C1::ScalarField::one()];
        let divisor_poly = U::from_coeffs(HostSlice::from_slice(&coeffs), 2);

        let (q_0, _) = &b.divide(&divisor_poly);
        
        let q_0 = Kzg::commit(pk, q_0).unwrap();
        let q_1 = Kzg::open(
            pk,
            &[witness.f.clone(), index.s.clone(), b, q],
            mu,
            separation_challenge,
        ).unwrap();

        Proof {
            sum,
            b_cm,
            q_cm,
            f_opening,
            s_opening,
            b_opening,
            q_opening,
            q_0,
            q_1,
        }
    }
    
    /// Verify the realtion if hold by reconstructing the equation.
    pub fn verify<Func>(
        index: &VerifierIndex<C1>,
        instance: &Instance<C1>,
        proof: &Proof<C1>,
        vk: &KzgVk<C1, C2, F>,
        relation: &Func,
    ) -> Result<(), Error>
    where
        Func: Fn(C1::ScalarField) -> C1::ScalarField,
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        let mut tr = Transcript::new_transcript(b"log-derivative");
        tr.send_v_index(index);
        tr.send_instance(instance);

        tr.send_blinders_sum(&proof.sum);
        let beta = tr.get_beta();
        let relation_at_beta = relation(beta);
        let b_0 =
            (relation_at_beta + proof.sum) * (C1::ScalarField::from_u32((N).try_into().unwrap()).inv());

        tr.send_b_and_q(&proof.b_cm, &proof.q_cm);
        let mu = tr.get_mu();

        tr.send_openings(
            &proof.f_opening,
            &proof.s_opening,
            &proof.b_opening,
            &proof.q_opening,
        );
        let separation_challenge = tr.get_separation_challenge();

        let sumcheck_relation = Kzg::verify(
            &[proof.b_cm],
            &[b_0],
            proof.q_0,
            C1::ScalarField::zero(),
            C1::ScalarField::one(),
            vk,
        );
        if sumcheck_relation.is_err() {
            return Err(Error::Sumcheck);
        }

        let openings_result = Kzg::verify(
            &[instance.f_cm, index.s_cm, proof.b_cm, proof.q_cm],
            &[
                proof.f_opening,
                proof.s_opening,
                proof.b_opening,
                proof.q_opening,
            ],
            proof.q_1,
            mu,
            separation_challenge,
            vk,
        );
        if openings_result.is_err() {
            return Err(Error::Openings);
        }
        
        // reconstruce the relation to checl if it hold.
        let formation_eq = {
            let zh_eval = mu.pow(N) - C1::ScalarField::one();
            proof.b_opening * (beta * proof.s_opening + proof.f_opening) - C1::ScalarField::one()
                == proof.q_opening * zh_eval
        };
        
        if !formation_eq {
            return Err(Error::WellFormation);
        }

        Ok(())
    }
}

#[cfg(test)]
mod log_derivative_tests {
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::Curve;
    use icicle_core::traits::FieldImpl;
    use std::marker::PhantomData;
    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, release_domain};
    use icicle_runtime::memory::HostSlice;
    use icicle_core::traits::Arithmetic;

    use crate::{
        kzg::{Kzg, PK, VK},
        utils::srs::unsafe_setup_from_tau,
    };

    use super::{
        structs::{Instance, Witness},
        Argument,
    };

    const N: usize = 16;
    const B: usize = 4;

    #[test]
    fn test_log_derivative() {
        let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

        let tau = Bn254ScalarField::from_u32(17u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData,};
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());

        let index_v = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::index_v(&pk);
        let index_p = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::index_p();

        // let's make f such that it has just one 1 and 14 zeros
        let mut f_evals = vec![Bn254ScalarField::zero(); N - B];
        f_evals[3] = Bn254ScalarField::one();

        let mut blinders: Vec<_> = (0..B).map(|i| Bn254ScalarField::from_u32((i + 10) as u32)).collect();
        let blinders_cloned = blinders.clone();
        f_evals.append(&mut blinders);
        
        let cfg = NTTConfig::<Bn254ScalarField>::default();
        let mut f_coeffs = vec![Bn254ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&f_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut f_coeffs),
        )
        .unwrap();

        let f = Bn254Poly::from_coeffs(HostSlice::from_slice(&f_coeffs), N);
        let f_cm = Kzg::commit(&pk, &f).unwrap();

        let instance = Instance::<Bn254CurveCfg> { f_cm };

        let witness = Witness::<Bn254CurveCfg, Bn254Poly> { f, e: PhantomData,};

        // RHS = 1/(beta + 1) + (N - B - 1)/(beta)
        let relation = |beta: Bn254ScalarField| {
            let beta_inv = beta.inv();
            let beta_plus_one_inv = (Bn254ScalarField::one() + beta).inv();
            let n_minus_one = Bn254ScalarField::from_u32((N - B - 1) as u32);

            beta_plus_one_inv + n_minus_one * beta_inv
        };

        let proof = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::prove(&index_p, &index_v, &instance, &witness, &pk);

        {
            let mut sum = Bn254ScalarField::zero();
            for i in 0..blinders_cloned.len() {
                sum = sum + blinders_cloned[i].inv();
            }
            assert_eq!(sum, proof.sum);
        }

        let result = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&index_v, &instance, &proof, &vk, &relation);

        release_domain::<Bn254ScalarField>().unwrap();

        assert!(result.is_ok());
    }
}