/*
    Given f(X) over domain of size N where b last evaluations are blinders
    we want to prove some log derivative relation "R(X)" about first N - b evaluations

    In other words: "prove that ∑1/(ß + f_i) = R(ß)"

    Both prover and verifier run indexer which computes s(X) such that
    s = [1, 1, ..., 0, 0, 0] (all 1 and then 0 last N - b evaluations)

    1. Prover sends ¥, which is claimed to be inverse blinders sum, ¥ = ∑1/fi for i in [N - b, N]
    2. Verifiers sends ß and computes R(ß)
    3. Prover sends b(X) and q(X)
    4. Verifier sends µ
    5. Prover sends b(0), b(µ), q(µ), f(µ), s(µ)
    6. Verifier sends separation challenge ¡
    7. Prover sends kzg proofs [πs]
    6. Verifier checks
        1. [πs]
        2. b(µ)(ß * s(µ) + f(µ)) - 1 = q(µ)zH(µ)
        3. b(0) = (R(ß) + ¥)/N
*/
use icicle_core::curve::Curve;
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, NTT, release_domain};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;
use crate::kzg::{Kzg, PK as KzgPk, VK as KzgVk};
use crate::utils::get_coeffs_of_poly;
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
    pub fn index_v(pk: &KzgPk<C1,C2,F>) -> VerifierIndex<C1> 
    where
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        //let domain = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());

        let mut zeros = vec![C1::ScalarField::zero(); B];
        let mut s_evals = vec![C1::ScalarField::one(); N - B];

        s_evals.append(&mut zeros);
         
        let cfg = NTTConfig::<C1::ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut coeffs = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&s_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut coeffs),
        )
        .unwrap();

        //release_domain::<C1::ScalarField>().unwrap();

        let s = U::from_coeffs(HostSlice::from_slice(&coeffs), N);
        let s_cm = Kzg::commit(pk, &s).unwrap();

        VerifierIndex { s_cm: s_cm.into() }
    }

    pub fn index_p() -> ProverIndex<C1, U> 
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        //let domain = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());
        let g = C1::ScalarField::from_u32(5u32);

        let mut twist = vec![C1::ScalarField::zero(); N];
        let mut pow = C1::ScalarField::one();
        for i in 0..N {
            twist[i] = pow;
            pow = pow * g;
        }

        let mut zeros = vec![C1::ScalarField::zero(); B];
        let mut s_evals = vec![C1::ScalarField::one(); N - B];

        s_evals.append(&mut zeros);
        
        let cfg = NTTConfig::<C1::ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut s_coeffs = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&s_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut s_coeffs),
        )
        .unwrap();

        let s_coeffs_clone = s_coeffs.clone();

        for i in 0..s_coeffs.len() {
            s_coeffs[i] = s_coeffs[i] * twist[i];
        }
        
        //domain
        let mut s_coset_evals = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&s_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut s_coset_evals),
        )
        .unwrap();
        
        //release_domain::<C1::ScalarField>().unwrap();

        ProverIndex {
            s: U::from_coeffs(HostSlice::from_slice(&s_coeffs_clone), N),
            s_coset_evals: s_coset_evals,
        }
    }

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
        //let domain = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());
        let g = C1::ScalarField::from_u32(5u32);

        let mut twist = vec![C1::ScalarField::zero(); N];
        let mut pow = C1::ScalarField::one();
        for i in 0..N {
            twist[i] = pow;
            pow = pow * g;
        }
        
        let mut f_coeffs = get_coeffs_of_poly(&witness.f);

        let cfg = NTTConfig::<C1::ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut f_evals = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&f_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut f_evals),
        )
        .unwrap();
        
        for i in 0..N {
            f_coeffs[i] = f_coeffs[i] * twist[i];
        }

        //domain
        let mut f_coset_evals = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&f_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut f_coset_evals),
        )
        .unwrap();

        let mut tr = Transcript::new(b"log-derivative");
        tr.send_v_index(v_index);
        tr.send_instance(instance);

        let mut blinders = f_evals[N - B..].to_vec();
        for i in 0..blinders.len() {
            blinders[i] = blinders[i].inv();
        }

        let mut gamma = C1::ScalarField::zero();
        for i in 0..blinders.len() {
            gamma = gamma + blinders[i];
        }

        tr.send_blinders_sum(&gamma);
        let beta = tr.get_beta();

        let mut b_evals = Vec::with_capacity(N - B);
        for i in 0..(N - B) {
            b_evals.push((f_evals[i] + beta).inv());
        }
        b_evals.append(&mut blinders);

        //domain
        let mut b_coeffs = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&b_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut b_coeffs),
        )
        .unwrap();

        let b = U::from_coeffs(HostSlice::from_slice(&b_coeffs), N);
        let b_cm = Kzg::commit(pk, &b).unwrap();
        
        for i in 0..N {
            b_coeffs[i] = b_coeffs[i] * twist[i];
        }

        //domain
        let mut b_coset_evals = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&b_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut b_coset_evals),
        )
        .unwrap();
        
        let x = C1::ScalarField::from_u32(5u32);
        let zh_coset_inv = (x.pow(N.try_into().unwrap()) - C1::ScalarField::one()).inv();

        let mut q_coset_evals = Vec::with_capacity(N);
        for i in 0..N {
            let b_i = b_coset_evals[i];
            let f_i = f_coset_evals[i];
            let s_i = index.s_coset_evals[i];

            let result = (b_i * (s_i * beta + f_i) - C1::ScalarField::one()) * zh_coset_inv;
            q_coset_evals.push(result);
        }
        
        //domain
        let mut q_coeffs = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&q_coset_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut q_coeffs),
        )
        .unwrap();

        for i in 0..N {
            q_coeffs[i] = q_coeffs[i] * twist[i].inv();
        }

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
        
        //release_domain::<C1::ScalarField>().unwrap();

        Proof {
            gamma,
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
        let mut tr = Transcript::new(b"log-derivative");
        tr.send_v_index(index);
        tr.send_instance(instance);

        tr.send_blinders_sum(&proof.gamma);
        let beta = tr.get_beta();
        let relation_at_beta = relation(beta);
        let b_0 =
            (relation_at_beta + proof.gamma) * (C1::ScalarField::from_u32((N).try_into().unwrap()).inv());

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

        //let domain = get_root_of_unity::<Bn254ScalarField>(N.try_into().unwrap());

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
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut f_coeffs = vec![Bn254ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&f_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut f_coeffs),
        )
        .unwrap();

        //release_domain::<Bn254ScalarField>().unwrap();

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

        /* */
        {
            let mut sum = Bn254ScalarField::zero();
            for i in 0..blinders_cloned.len() {
                sum = sum + blinders_cloned[i].inv();
            }
            assert_eq!(sum, proof.gamma);
        }

        let result = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&index_v, &instance, &proof, &vk, &relation);

        release_domain::<Bn254ScalarField>().unwrap();

        assert!(result.is_ok());
    }
}