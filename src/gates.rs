//! Gates argument over a 2N-domain with KZG commitments and NTT-based evaluation.
//!
//! This module builds and verifies a 4-gate arithmetic circuit used to tie together
//! public/hidden bids with randomizers. All constraints are aggregated into a single
//! quotient polynomial `q(X)` on an extended coset and committed with KZG.
//!
//! # Gates
//! - **G1 (Inverse):** `q_price(X) · (r(X)·r_inv(X) − 1) = 0`
//! - **G2 (Link):** `q_price(X) · (g(X) − f(X) − bid(X)·r(X)) = 0`
//! - **G3 (Shift):** `q_price(X) · (diff(X) − bid(X) + bid(ωX)) = 0`
//! - **G4 (Boundary):** `L_p(X) · bid(X) = 0`
//!
//! Here `q_price` selects the active price rows.
//!
//! degree of quotient will be n - 1 + n - 1 + n - 1 - n = 3n - 3 - n = 2n - 3
use icicle_core::curve::Curve;
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, ntt_inplace, NTT, release_domain};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;
use crate::utils::{my_ntt, get_device_is_cpu_or_gpu};

use crate::{
    kzg::{Kzg, PK as KzgPk, VK as KzgVk},
    utils::{get_coeffs_of_poly, evaluate_vanishing_over_extended_coset},
};

use self::{
    structs::{Error, Oracle, Proof, ProverIndex, VerifierIndex, Witness},
    tr::Transcript,
};

pub mod structs;
mod tr;

impl<'a, F: FieldImpl> Oracle<'a, F> {
    pub fn query(&self, i: usize, rotation: usize, extension: usize) -> F {
        self.0[(i + rotation * extension) % self.0.len()]
    }
}

pub struct GatesArgument<const N: usize, const P: usize, C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    _e: PhantomData<(C1, C2, F)>
}

impl<const N: usize, const P: usize, C1, C2, F> GatesArgument<N, P, C1, C2, F>
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    /// A helper function to constructe the selector vector
    fn make_q_price<U>() -> U
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
     {
        let cpu_or_gpu = get_device_is_cpu_or_gpu();

        let mut q_price_evals = vec![C1::ScalarField::one(); P];
        let mut zeros = vec![C1::ScalarField::zero(); N - P];
        q_price_evals.append(&mut zeros);
        
        let coeffs = my_ntt::<C1>(&q_price_evals, C1::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);

        U::from_coeffs(HostSlice::from_slice(&coeffs), N) 
    }
    
    /// A helper function for verifier to constructe the selector vector
    pub fn verifier_index<U>(pk: &KzgPk<C1,C2,F>) -> VerifierIndex<C1>
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        let q_price = Self::make_q_price::<U>();
        let q_price_cm = Kzg::commit(pk, &q_price).unwrap();
        VerifierIndex { q_price_cm }
    }
   
    /// A helper function for the prover to construct get the vectors for proving.
    pub fn prover_index<U>() -> ProverIndex<C1, U> 
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic,
     {  
        let cpu_or_gpu = get_device_is_cpu_or_gpu();
        let g = C1::ScalarField::from_u32(5u32);
        
        let q_price: U = Self::make_q_price();
        
        let mut q_price_coeffs = get_coeffs_of_poly(&q_price);
        q_price_coeffs.resize(2 * N, C1::ScalarField::zero());
        
        //domain_2n
        let q_price_coset_evals = my_ntt::<C1>(&q_price_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        
        let mut l_p_evals = vec![C1::ScalarField::zero(); N];
        l_p_evals[P] = C1::ScalarField::one();
        
        //domain_n
        let mut l_p = my_ntt::<C1>(&l_p_evals, C1::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        l_p.resize(2 * N, C1::ScalarField::zero());
        
        //domain_2n
        let mut l_p_coset_evals = my_ntt::<C1>(&l_p, g, NTTDir::kForward, cpu_or_gpu);
        
        ProverIndex {
            q_price,
            q_price_coset_evals,
            l_p_coset_evals,
        }
    }
    
    // The function is used to construct the 4 gates and compute the quotient poly of the sum of the 4 gates.
    pub fn prove<U>(
        witness: &Witness<C1, U>,
        v_index: &VerifierIndex<C1>,
        index: &ProverIndex<C1, U>,
        pk: &KzgPk<C1, C2, F>,
    ) -> Proof<C1>
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        let cpu_or_gpu = get_device_is_cpu_or_gpu();
        let k = 2;
        
        let domain_n = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());

        let mut tr = Transcript::<C1>::new_transcript(b"gates-transcript");
        tr.send_index(v_index);

        let bid_cm = Kzg::commit(pk, &witness.bid).unwrap();
        let random_x_cm = Kzg::commit(pk, &witness.random_x).unwrap();
        let random_r_cm = Kzg::commit(pk, &witness.random_r).unwrap();
        let random_r_inv_cm = Kzg::commit(pk, &witness.random_r_inv).unwrap();
        let diff_f_cm = Kzg::commit(pk, &witness.diff_f).unwrap();
        let hidden_bid_cm = Kzg::commit(pk, &witness.hidden_bid).unwrap();
        
        tr.send_oracle_commitments(&bid_cm, &random_x_cm, &random_r_cm, &random_r_inv_cm, &diff_f_cm, &hidden_bid_cm);
        let alpha = tr.get_quotient_challenge();
        let alpha_pows: Vec<C1::ScalarField> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * alpha)).take(4).collect();

        let g = C1::ScalarField::from_u32(5u32);
        
        let mut bid_coeffs = get_coeffs_of_poly(&witness.bid);
        bid_coeffs.resize(k * N, C1::ScalarField::zero());

        //domain_kn
        let bid_coset_evals = my_ntt::<C1>(&bid_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        let bid_coset_evals = Oracle(&bid_coset_evals);
        
        let mut random_x_coeffs = get_coeffs_of_poly(&witness.random_x);
        random_x_coeffs.resize(k * N, C1::ScalarField::zero());

        //domain_kn
        let random_x_coset_evals = my_ntt::<C1>(&random_x_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        let random_x_coset_evals = Oracle(&random_x_coset_evals);

        let mut random_r_coeffs = get_coeffs_of_poly(&witness.random_r);
        random_r_coeffs.resize(k * N, C1::ScalarField::zero());

        //domain_kn
        let random_r_coset_evals = my_ntt::<C1>(&random_r_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        let random_r_coset_evals = Oracle(&random_r_coset_evals);

        let mut random_r_inv_coeffs = get_coeffs_of_poly(&witness.random_r_inv);
        random_r_inv_coeffs.resize(k * N, C1::ScalarField::zero());
        
        //domain_kn
        let random_r_inv_coset_evals = my_ntt::<C1>(&random_r_inv_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        let random_r_inv_coset_evals = Oracle(&random_r_inv_coset_evals);

        let mut diff_f_coeffs = get_coeffs_of_poly(&witness.diff_f);
        diff_f_coeffs.resize(k * N, C1::ScalarField::zero());

        //domain_kn
        let diff_f_coset_evals = my_ntt::<C1>(&diff_f_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        let diff_f_coset_evals = Oracle(&diff_f_coset_evals);

        let mut hidden_bid_coeffs = get_coeffs_of_poly(&witness.hidden_bid);
        hidden_bid_coeffs.resize(k * N, C1::ScalarField::zero());

        //domain_kn
        let hidden_bid_coset_evals = my_ntt::<C1>(&hidden_bid_coeffs, g, NTTDir::kForward, cpu_or_gpu);
        let hidden_bid_coset_evals = Oracle(&hidden_bid_coset_evals);

        let q_price_coset_evals = Oracle(&index.q_price_coset_evals);
        let l_p_coset_evals = Oracle(&index.l_p_coset_evals);
        
        let modulus_zh_coset_evals =
            evaluate_vanishing_over_extended_coset::<C1>(N, k);
        
        let mut modulus_zh_coset_evals_inv = Vec::new();
        for i in &modulus_zh_coset_evals {
            modulus_zh_coset_evals_inv.push(i.inv());
        }

        let mut q_coset_evals = vec![C1::ScalarField::zero(); k * N];
        let one = C1::ScalarField::one();
        
        // computes the qutotient poly for the four gates
        let compute_q_for = |i: usize| -> <C1 as Curve>::ScalarField {
            let q_price_i = q_price_coset_evals.query(i, 0, k);
            let bid_i = bid_coset_evals.query(i, 0, k);
            let bid_i_next = bid_coset_evals.query(i, 1, k); // value at w·X
            let random_r_i = random_r_coset_evals.query(i, 0, k);
            let random_r_inv_i = random_r_inv_coset_evals.query(i, 0, k);
            let random_x_i = random_x_coset_evals.query(i, 0, k);
            let diff_f_i  = diff_f_coset_evals.query(i, 0, k);
            let hidden_bid_i = hidden_bid_coset_evals.query(i, 0, k);
            let l_p_i = l_p_coset_evals.query(i, 0, k);

            // Gate 1: q_price * (r * r_inv - 1)
            let g1 = alpha_pows[0] * q_price_i * (random_r_i * random_r_inv_i - one);

            // Gate 2: q_price * (g - f - bid * r)
            let g2 = alpha_pows[1] * q_price_i * (hidden_bid_i - random_x_i - bid_i * random_r_i);

            // Gate 3: q_price * (diff - bid + bid(w·X))
            let g3 = alpha_pows[2] * q_price_i * (diff_f_i - bid_i + bid_i_next);

            // Gate 4: L_p * bid  (activates only at position p)
            let g4 = alpha_pows[3] * l_p_i * bid_i;

            // Rescale by the vanishing polynomial on the extended domain.
            let zh_inv_i = modulus_zh_coset_evals_inv[i % k];

            (g1 + g2 + g3 + g4) * zh_inv_i
        };

        // Fill q_coset_evals using the closure above.
        q_coset_evals
            .iter_mut()
            .enumerate()
            .for_each(|(i, q)| *q = compute_q_for(i));
        
        //domain_kn
        let q = my_ntt::<C1>(&q_coset_evals, g, NTTDir::kInverse, cpu_or_gpu);

        // hardcode to 2 chunks
        let q_chunk_0 = U::from_coeffs(HostSlice::from_slice(&q[..N]), q[..N].len());
        let q_chunk_1 = U::from_coeffs(HostSlice::from_slice(&q[N..]), q[N..].len());

        let q_chunk_0_cm = Kzg::commit(pk, &q_chunk_0).unwrap();
        let q_chunk_1_cm = Kzg::commit(pk, &q_chunk_1).unwrap();

        tr.send_q_chunks(&q_chunk_0_cm, &q_chunk_1_cm);

        // open everything
        let gamma = tr.get_evaluation_challenge();
        
        let q_price_opening = index.q_price.eval(&gamma);
        let bid_opening = witness.bid.eval(&gamma);
        let bid_shift_opening = witness.bid.eval(&(gamma * domain_n));
        let random_x_opening = witness.random_x.eval(&gamma);
        let random_r_opening = witness.random_r.eval(&gamma);
        let random_r_inv_opening = witness.random_r_inv.eval(&gamma);
        let diff_f_opening = witness.diff_f.eval(&gamma);
        let hidden_bid_opening = witness.hidden_bid.eval(&gamma);
        let q_chunk_0_opening = q_chunk_0.eval(&gamma);
        let q_chunk_1_opening = q_chunk_1.eval(&gamma);

        tr.send_oracle_openings(
            &q_price_opening,
            &bid_opening,
            &bid_shift_opening,
            &random_x_opening,
            &random_r_opening,
            &random_r_inv_opening,
            &diff_f_opening,
            &hidden_bid_opening,
            &q_chunk_0_opening,
            &q_chunk_1_opening,
        );

        let separation_challenge = tr.get_separation_challenge();
        
        let w_0 = Kzg::open(
            pk,
            &[
                index.q_price.clone(),
                witness.bid.clone(),
                witness.random_x.clone(),
                witness.random_r.clone(),
                witness.random_r_inv.clone(),
                witness.diff_f.clone(),
                witness.hidden_bid.clone(),
                q_chunk_0.clone(),
                q_chunk_1.clone(),
            ],
            gamma,
            separation_challenge,
        ).unwrap();

        let w_1 = Kzg::open(pk, &[witness.bid.clone()], gamma * domain_n, one).unwrap();

        Proof {
            bid_cm,
            random_r_cm,
            random_r_inv_cm,
            random_x_cm,
            diff_f_cm,
            hidden_bid_cm,
            q_price_opening,
            bid_opening,
            bid_shift_opening,
            random_x_opening,
            random_r_opening,
            random_r_inv_opening,
            diff_f_opening,
            hidden_bid_opening,
            q_chunk_0_cm,
            q_chunk_1_cm,
            q_chunk_0_opening,
            q_chunk_1_opening,
            w_0,
            w_1,
        }
    }
    
    // Check the 4 gates relation hold by reconstruct it.
    pub fn verify(
        index: &VerifierIndex<C1>,
        proof: &Proof<C1>,
        vk: &KzgVk<C1, C2, F>,
    ) -> Result<(), Error> 
    where
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic
    {
        let domain = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());
        let mut tr = Transcript::<C1>::new_transcript(b"gates-transcript");
        tr.send_index(index);

        tr.send_oracle_commitments(
            &proof.bid_cm,
            &proof.random_x_cm,
            &proof.random_r_cm,
            &proof.random_r_inv_cm,
            &proof.diff_f_cm,
            &proof.hidden_bid_cm,
        );
        let alpha = tr.get_quotient_challenge();
        let alpha_pows: Vec<C1::ScalarField> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * alpha)).take(4).collect();

        tr.send_q_chunks(&proof.q_chunk_0_cm, &proof.q_chunk_1_cm);
        let gamma = tr.get_evaluation_challenge();
        
        tr.send_oracle_openings(
            &proof.q_price_opening,
            &proof.bid_opening,
            &proof.bid_shift_opening,
            &proof.random_x_opening,
            &proof.random_r_opening,
            &proof.random_r_inv_opening,
            &proof.diff_f_opening,
            &proof.hidden_bid_opening,
            &proof.q_chunk_0_opening,
            &proof.q_chunk_1_opening,
        );

        let separation_challenge = tr.get_separation_challenge();

        let res_gamma = Kzg::verify(
            &[
                index.q_price_cm,
                proof.bid_cm,
                proof.random_x_cm,
                proof.random_r_cm,
                proof.random_r_inv_cm,
                proof.diff_f_cm,
                proof.hidden_bid_cm,
                proof.q_chunk_0_cm,
                proof.q_chunk_1_cm,
            ],
            &[
                proof.q_price_opening,
                proof.bid_opening,
                proof.random_x_opening,
                proof.random_r_opening,
                proof.random_r_inv_opening,
                proof.diff_f_opening,
                proof.hidden_bid_opening,
                proof.q_chunk_0_opening,
                proof.q_chunk_1_opening,
            ],
            proof.w_0,
            gamma,
            separation_challenge,
            vk,
        );

        if !res_gamma.is_ok() {
            return Err(Error::Opening);
        }

        let res_gamma_sh = Kzg::verify(
            &[proof.bid_cm],
            &[proof.bid_shift_opening],
            proof.w_1,
            gamma * domain,
            C1::ScalarField::one(),
            vk,
        );

        if !res_gamma_sh.is_ok() {
            return Err(Error::ShiftedOpening);
        }

        let zh_at_gamma = gamma.pow(N) - C1::ScalarField::one();
        let l_p_next_at_gamma = {
            let n = C1::ScalarField::from_u32(N as u32);
            let n_inv = n.inv();
            let w_p = domain.pow((P).try_into().unwrap());

            let x_minus_w_p_inv = (gamma - w_p).inv();

            w_p * n_inv * zh_at_gamma * x_minus_w_p_inv
        };

        let gamma_pow_n = gamma.pow(N);
        
        // check the 4 gates relation
        let lhs = {
            let g1 = alpha_pows[0]
                * proof.q_price_opening
                * (proof.random_r_opening * proof.random_r_inv_opening - C1::ScalarField::one());

            let g2 = alpha_pows[1]
                * proof.q_price_opening
                * (proof.hidden_bid_opening - proof.random_x_opening - proof.bid_opening * proof.random_r_opening);

            let g3 = alpha_pows[2]
                * proof.q_price_opening
                * (proof.diff_f_opening - proof.bid_opening + proof.bid_shift_opening);

            let g4 = alpha_pows[3] * l_p_next_at_gamma * proof.bid_opening;

            g1 + g2 + g3 + g4
        };

        let rhs =
            { (proof.q_chunk_0_opening + gamma_pow_n * proof.q_chunk_1_opening) * zh_at_gamma };

        if lhs != rhs {
            return Err(Error::RelationCheck);
        }

        Ok(())
    }
}

#[cfg(test)]
mod gates_test {
    use std::ops::Mul;
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::Curve;
    use icicle_core::traits::FieldImpl;
    use rand_chacha::ChaCha20Rng;
    use std::marker::PhantomData;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;

    use crate::bid_encoder::BidEncoder;
    use crate::{
        kzg::{PK, VK},
        utils::srs::unsafe_setup_from_tau,
    };

    use super::GatesArgument;

    use icicle_core::ntt::{NTTInitDomainConfig, initialize_domain, release_domain, get_root_of_unity};

    const P: usize = 10;
    const N: usize = 16;

    const SEED: [u8; 32] = [
        1, 0, 52, 0, 0, 0, 0, 0, 1, 0, 10, 0, 22, 32, 0, 0, 2, 0, 55, 49, 0, 11, 0, 0, 3, 0, 0, 0,
        0, 0, 2, 92,
    ];

    #[test]
    fn test_gates() {
        type Poly = Bn254Poly;

        let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        
        let tau = Bn254ScalarField::from_u32(17u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator().mul(tau);

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());
        
        let v_index = GatesArgument::<N, P, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::verifier_index::<Poly>(&pk);
        let p_index = GatesArgument::<N, P, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::prover_index::<Poly>();
        
        let bid = 9;
        let enc = BidEncoder::<P, N, Bn254CurveCfg>::encode::<ChaCha20Rng>(bid, SEED);
        
        let witness= enc.to_gate_witness::<ChaCha20Rng, Poly>(SEED);
        
        let proof = GatesArgument::<N, P, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::prove(&witness, &v_index, &p_index, &pk);
        let result = GatesArgument::<N, P, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::verify(&v_index, &proof, &vk);
        
        release_domain::<Bn254ScalarField>().unwrap();
        
        assert!(result.is_ok());
    }
}

#[cfg(test)]
mod ntt_test {
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;
    use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, release_domain};
    use icicle_runtime::memory::HostSlice;
    use icicle_core::traits::Arithmetic;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;
    
    const N: usize = 16;

    #[test]
    #[ignore]
    fn test_my_ntt() {
        let domain_2n = get_root_of_unity::<Bn254ScalarField>((2 * N).try_into().unwrap());
        initialize_domain(domain_2n, &NTTInitDomainConfig::default()).unwrap();
        
        let g = Bn254ScalarField::from_u32(5u32);
        let cfg = NTTConfig::<Bn254ScalarField>::default();
        let mut cfg_2n = NTTConfig::<Bn254ScalarField>::default();
        cfg_2n.coset_gen = g;
        
        let mut twist = vec![Bn254ScalarField::zero(); 2 * N];
        let mut pow = Bn254ScalarField::one();
        for i in 0..(2 * N) {
            twist[i] = pow;
            pow = pow * g;
        }
        
        let target = vec![ScalarCfg::generate_random(1)[0] ; N];
        let mut target_copy = target.clone();
        let mut target_copy_2 = target.clone();

        target_copy.resize(2 * N, Bn254ScalarField::zero());
        target_copy_2.resize(2 * N, Bn254ScalarField::zero());

        for i in 0..(2 * N) {
            target_copy_2[i] = twist[i] * target_copy_2[i];
        }

        for i in 0..(2 * N) {
            target_copy_2[i] = twist[i].inv() * target_copy_2[i];
        }

        assert_eq!(target_copy_2, target_copy);

        for i in 0..(2 * N) {
            target_copy_2[i] = twist[i] * target_copy_2[i];
        }

        let mut q = vec![Bn254ScalarField::zero(); target.len()];
        //domain_2n
        ntt(
            HostSlice::from_slice(&target),
            NTTDir::kForward,
            &cfg_2n,
            HostSlice::from_mut_slice(&mut q),
        )
        .unwrap();
        
        let mut q_copy = vec![Bn254ScalarField::zero(); target_copy.len()];
        //domain_2n
        ntt(
            HostSlice::from_slice(&target_copy),
            NTTDir::kForward,
            &cfg_2n,
            HostSlice::from_mut_slice(&mut q_copy),
        )
        .unwrap();
        
        let mut q_g = vec![Bn254ScalarField::zero(); target_copy.len()];
        //domain_2n
        ntt(
            HostSlice::from_slice(&target_copy_2),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut q_g),
        )
        .unwrap();
        assert_eq!(q_g, q_copy);

        let mut q_copy_inv = vec![Bn254ScalarField::zero(); target_copy.len()];
        //domain_2n
        ntt(
            HostSlice::from_slice(&q_copy),
            NTTDir::kInverse,
            &cfg_2n,
            HostSlice::from_mut_slice(&mut q_copy_inv),
        )
        .unwrap();

        assert_eq!(q_copy_inv, target_copy);

        let mut q_g_inv = vec![Bn254ScalarField::zero(); q_g.len()];
        //domain_2n
        ntt(
            HostSlice::from_slice(&q_g),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut q_g_inv),
        )
        .unwrap();
        
        for i in 0..(2 * N) {
            q_g_inv[i] = twist[i].inv() * q_g_inv[i];
        }

        release_domain::<Bn254ScalarField>().unwrap();

        assert_eq!(q_g_inv, target_copy);
    }
}