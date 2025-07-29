/*
    All gates exist over q_price(X) which is selector that is 1 on |p| - price range entires
    gate1: q_price(X)(r(X) * r_inv(X) - 1) = 0 mod zH(X)
    gate2: q_price(X)(g(X) - f(X) - bid(X)*r(X)) = 0 mod zH(X)
    gate3: q_price(X)(diff(X) - bid(X) + bid(wX)) = 0 mod zH(X)
    gate4: L_p(X)bid(X) = 0 mod zH(X)

    degree of quotient will be n - 1 + n - 1 + n - 1 - n = 3n - 3 - n = 2n - 3, so we can work with subgroup of 2n
*/
use icicle_core::curve::Curve;
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, NTT, release_domain};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;

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
    fn make_q_price<U>() -> U
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
     {
        let domain = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());

        let mut q_price_evals = vec![C1::ScalarField::one(); P];
        let mut zeros = vec![C1::ScalarField::zero(); N - P];
        q_price_evals.append(&mut zeros);
        
        let cfg = NTTConfig::<C1::ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut coeffs = vec![C1::ScalarField::zero(); N];
        ntt(
            HostSlice::from_slice(&q_price_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut coeffs),
        )
        .unwrap();
        
        //release_domain::<C1::ScalarField>().unwrap();

        U::from_coeffs(HostSlice::from_slice(&coeffs), N) 
    }

    pub fn verifier_index<U>(pk: &KzgPk<C1,C2,F>) -> VerifierIndex<C1>
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        let q_price = Self::make_q_price::<U>();
        let q_price_cm = Kzg::commit(pk, &q_price);
        VerifierIndex { q_price_cm }
    }

    pub fn prover_index<U>() -> ProverIndex<C1, U> 
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic,
     {  
        let g = C1::ScalarField::from_u32(5u32);
        let mut twist = vec![C1::ScalarField::zero(); 2 * N];
        let mut pow = C1::ScalarField::one();
        
        for i in 0..(2 * N) {
            twist[i] = pow;
            pow = pow * g;
        }
        
        let q_price: U = Self::make_q_price();
        
        let mut q_price_coeffs = get_coeffs_of_poly(&q_price);
        q_price_coeffs.resize(2 * N, C1::ScalarField::zero());
        for i in 0..(2 * N) {
            q_price_coeffs[i] = q_price_coeffs[i] * twist[i];
        }

        let cfg = NTTConfig::<C1::ScalarField>::default();
        let domain_2n = get_root_of_unity::<C1::ScalarField>((2 * N).try_into().unwrap());
        //initialize_domain(domain_2n, &NTTInitDomainConfig::default()).unwrap();
        
        let mut q_price_coset_evals = vec![C1::ScalarField::zero(); 2 * N];
        //domain_2n
        ntt(
            HostSlice::from_slice(&q_price_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut q_price_coset_evals),
        )
        .unwrap();
        
        let mut l_p_evals = vec![C1::ScalarField::zero(); N];
        l_p_evals[P] = C1::ScalarField::one();

        let mut l_p = vec![C1::ScalarField::zero(); N];
        
        //domain_n
        ntt(
            HostSlice::from_slice(&l_p_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut l_p),
        )
        .unwrap();
        
        let mut l_p_coset_evals = vec![C1::ScalarField::zero(); 2 * N];
        l_p.resize(2 * N, C1::ScalarField::zero());
        for i in 0..(2 * N) {
            l_p[i] = l_p[i] * twist[i];
        }

        //domain_2n
        ntt(
            HostSlice::from_slice(&l_p),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut l_p_coset_evals),
        )
        .unwrap();
        
        //release_domain::<C1::ScalarField>().unwrap();
        
        ProverIndex {
            q_price,
            q_price_coset_evals,
            l_p_coset_evals,
        }
    }

    pub fn prove<U>(
        witness: &Witness<C1, U>,
        v_index: &VerifierIndex<C1>, // just to hash
        index: &ProverIndex<C1, U>,
        pk: &KzgPk<C1, C2, F>,
    ) -> Proof<C1>
    where
        U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: FieldImpl,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic,
    {
        let k = 2;
        
        let domain_n = get_root_of_unity::<C1::ScalarField>(N.try_into().unwrap());

        let mut tr = Transcript::<C1>::new(b"gates-transcript");
        tr.send_index(v_index);

        let bid_cm = Kzg::commit(pk, &witness.bid);
        let f_cm = Kzg::commit(pk, &witness.f);
        let r_cm = Kzg::commit(pk, &witness.r);
        let r_inv_cm = Kzg::commit(pk, &witness.r_inv);
        let diff_cm = Kzg::commit(pk, &witness.diff);
        let g_cm = Kzg::commit(pk, &witness.g);

        tr.send_oracle_commitments(&bid_cm, &f_cm, &r_cm, &r_inv_cm, &diff_cm, &g_cm);
        let alpha = tr.get_quotient_challenge();
        let alpha_pows: Vec<C1::ScalarField> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * alpha)).take(4).collect();

        let g = C1::ScalarField::from_u32(5u32);
        let mut twist = vec![C1::ScalarField::zero(); 2 * N];
        let mut pow = C1::ScalarField::one();
        for i in 0..(2 * N) {
            twist[i] = pow;
            pow = pow * g;
        }
        
        let mut bid_coeffs = get_coeffs_of_poly(&witness.bid);
        bid_coeffs.resize(k * N, C1::ScalarField::zero());
        let mut bid_coset_evals = vec![C1::ScalarField::zero(); k * N];

        for i in 0..(k * N) {
            bid_coeffs[i] = bid_coeffs[i] * twist[i];
        }
        
        let cfg = NTTConfig::<C1::ScalarField>::default();
        let domain_kn = get_root_of_unity::<C1::ScalarField>((k * N).try_into().unwrap());
        //initialize_domain(domain_kn, &NTTInitDomainConfig::default()).unwrap();
        //domain_kn
        ntt(
            HostSlice::from_slice(&bid_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut bid_coset_evals),
        )
        .unwrap();
        
        let bid_coset_evals = Oracle(&bid_coset_evals);
        
        let mut f_coeffs = get_coeffs_of_poly(&witness.f);
        f_coeffs.resize(k * N, C1::ScalarField::zero());
        let mut f_coset_evals = vec![C1::ScalarField::zero(); k * N];

        for i in 0..(k * N) {
            f_coeffs[i] = f_coeffs[i] * twist[i];
        }

        //domain_kn
        ntt(
            HostSlice::from_slice(&f_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut f_coset_evals),
        )
        .unwrap();
        let f_coset_evals = Oracle(&f_coset_evals);

        let mut r_coeffs = get_coeffs_of_poly(&witness.r);
        r_coeffs.resize(k * N, C1::ScalarField::zero());
        let mut r_coset_evals = vec![C1::ScalarField::zero(); k * N];

        for i in 0..(k * N) {
            r_coeffs[i] = r_coeffs[i] * twist[i];
        }

        //domain_kn
        ntt(
            HostSlice::from_slice(&r_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut r_coset_evals),
        )
        .unwrap();
        let r_coset_evals = Oracle(&r_coset_evals);

        let mut r_inv_coeffs = get_coeffs_of_poly(&witness.r_inv);
        r_inv_coeffs.resize(k * N, C1::ScalarField::zero());
        let mut r_inv_coset_evals = vec![C1::ScalarField::zero(); k * N];

        for i in 0..(k * N) {
            r_inv_coeffs[i] = r_inv_coeffs[i] * twist[i];
        }

        //domain_kn
        ntt(
            HostSlice::from_slice(&r_inv_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut r_inv_coset_evals),
        )
        .unwrap();
        let r_inv_coset_evals = Oracle(&r_inv_coset_evals);

        let mut diff_coeffs = get_coeffs_of_poly(&witness.diff);
        diff_coeffs.resize(k * N, C1::ScalarField::zero());
        let mut diff_coset_evals = vec![C1::ScalarField::zero(); k * N];

        for i in 0..(k * N) {
            diff_coeffs[i] = diff_coeffs[i] * twist[i];
        }

        //domain_kn
        ntt(
            HostSlice::from_slice(&diff_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut diff_coset_evals),
        )
        .unwrap();
        let diff_coset_evals = Oracle(&diff_coset_evals);

        let mut g_coeffs = get_coeffs_of_poly(&witness.g);
        g_coeffs.resize(k * N, C1::ScalarField::zero());
        let mut g_coset_evals = vec![C1::ScalarField::zero(); k * N];

        for i in 0..(k * N) {
            g_coeffs[i] = g_coeffs[i] * twist[i];
        }

        //domain_kn
        ntt(
            HostSlice::from_slice(&g_coeffs),
            NTTDir::kForward,
            &cfg,
            HostSlice::from_mut_slice(&mut g_coset_evals),
        )
        .unwrap();
        let g_coset_evals = Oracle(&g_coset_evals);

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

        for i in 0..(k * N) {
            let q_price_i = q_price_coset_evals.query(i, 0, k);
            let bid_i = bid_coset_evals.query(i, 0, k);
            let bid_i_next = bid_coset_evals.query(i, 1, k);
            let r_i = r_coset_evals.query(i, 0, k);
            let r_inv_i = r_inv_coset_evals.query(i, 0, k);
            let f_i = f_coset_evals.query(i, 0, k);
            let diff_i = diff_coset_evals.query(i, 0, k);
            let g_i = g_coset_evals.query(i, 0, k);
            let l_p_i = l_p_coset_evals.query(i, 0, k);

            // gate1
            q_coset_evals[i] = alpha_pows[0] * q_price_i * (r_i * r_inv_i - one);

            // gate2
            q_coset_evals[i] = q_coset_evals[i] + (alpha_pows[1] * q_price_i * (g_i - f_i - bid_i * r_i));

            // gate3
            q_coset_evals[i] = q_coset_evals[i] + (alpha_pows[2] * q_price_i * (diff_i - bid_i + bid_i_next));

            // gate4
            q_coset_evals[i] = q_coset_evals[i] + (alpha_pows[3] * l_p_i * bid_i);

            // rescale by zh_inv
            let zh_inv_i = modulus_zh_coset_evals_inv[i % k];
            q_coset_evals[i] = q_coset_evals[i] * zh_inv_i
        }

        let mut q = vec![C1::ScalarField::zero(); q_coset_evals.len()];
        //domain_kn
        ntt(
            HostSlice::from_slice(&q_coset_evals),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut q),
        )
        .unwrap();
        
        //release_domain::<C1::ScalarField>().unwrap();

        for i in 0..q.len() {
            q[i] = q[i] * (twist[i].inv());
        }
 
        // hardcoded to 2 chunks
        let q_chunk_0 = U::from_coeffs(HostSlice::from_slice(&q[..N]), q[..N].len());
        let q_chunk_1 = U::from_coeffs(HostSlice::from_slice(&q[N..]), q[N..].len());

        let q_chunk_0_cm = Kzg::commit(pk, &q_chunk_0);
        let q_chunk_1_cm = Kzg::commit(pk, &q_chunk_1);

        tr.send_q_chunks(&q_chunk_0_cm, &q_chunk_1_cm);

        // open everything
        let gamma = tr.get_evaluation_challenge();
        
        let q_price_opening = index.q_price.eval(&gamma);
        let bid_opening = witness.bid.eval(&gamma);
        let bid_shift_opening = witness.bid.eval(&(gamma * domain_n));
        let f_opening = witness.f.eval(&gamma);
        let r_opening = witness.r.eval(&gamma);
        let r_inv_opening = witness.r_inv.eval(&gamma);
        let diff_opening = witness.diff.eval(&gamma);
        let g_opening = witness.g.eval(&gamma);
        let q_chunk_0_opening = q_chunk_0.eval(&gamma);
        let q_chunk_1_opening = q_chunk_1.eval(&gamma);

        tr.send_oracle_openings(
            &q_price_opening,
            &bid_opening,
            &bid_shift_opening,
            &f_opening,
            &r_opening,
            &r_inv_opening,
            &diff_opening,
            &g_opening,
            &q_chunk_0_opening,
            &q_chunk_1_opening,
        );

        let separation_challenge = tr.get_separation_challenge();

        let w_0 = Kzg::open(
            pk,
            &[
                index.q_price.clone(),
                witness.bid.clone(),
                witness.f.clone(),
                witness.r.clone(),
                witness.r_inv.clone(),
                witness.diff.clone(),
                witness.g.clone(),
                q_chunk_0.clone(),
                q_chunk_1.clone(),
            ],
            gamma,
            separation_challenge,
        );

        let w_1 = Kzg::open(pk, &[witness.bid.clone()], gamma * domain_n, one);

        Proof {
            bid_cm,
            r_cm,
            r_inv_cm,
            f_cm,
            diff_cm,
            g_cm,
            q_price_opening,
            bid_opening,
            bid_shift_opening,
            f_opening,
            r_opening,
            r_inv_opening,
            diff_opening,
            g_opening,
            q_chunk_0_cm,
            q_chunk_1_cm,
            q_chunk_0_opening,
            q_chunk_1_opening,
            w_0,
            w_1,
        }
    }

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
        let mut tr = Transcript::<C1>::new(b"gates-transcript");
        tr.send_index(index);

        tr.send_oracle_commitments(
            &proof.bid_cm,
            &proof.f_cm,
            &proof.r_cm,
            &proof.r_inv_cm,
            &proof.diff_cm,
            &proof.g_cm,
        );
        let alpha = tr.get_quotient_challenge();
        let alpha_pows: Vec<C1::ScalarField> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * alpha)).take(4).collect();

        tr.send_q_chunks(&proof.q_chunk_0_cm, &proof.q_chunk_1_cm);
        let gamma = tr.get_evaluation_challenge();
        
        tr.send_oracle_openings(
            &proof.q_price_opening,
            &proof.bid_opening,
            &proof.bid_shift_opening,
            &proof.f_opening,
            &proof.r_opening,
            &proof.r_inv_opening,
            &proof.diff_opening,
            &proof.g_opening,
            &proof.q_chunk_0_opening,
            &proof.q_chunk_1_opening,
        );

        let separation_challenge = tr.get_separation_challenge();

        let res_gamma = Kzg::verify(
            &[
                index.q_price_cm,
                proof.bid_cm,
                proof.f_cm,
                proof.r_cm,
                proof.r_inv_cm,
                proof.diff_cm,
                proof.g_cm,
                proof.q_chunk_0_cm,
                proof.q_chunk_1_cm,
            ],
            &[
                proof.q_price_opening,
                proof.bid_opening,
                proof.f_opening,
                proof.r_opening,
                proof.r_inv_opening,
                proof.diff_opening,
                proof.g_opening,
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

        let lhs = {
            let g1 = alpha_pows[0]
                * proof.q_price_opening
                * (proof.r_opening * proof.r_inv_opening - C1::ScalarField::one());

            let g2 = alpha_pows[1]
                * proof.q_price_opening
                * (proof.g_opening - proof.f_opening - proof.bid_opening * proof.r_opening);

            let g3 = alpha_pows[2]
                * proof.q_price_opening
                * (proof.diff_opening - proof.bid_opening + proof.bid_shift_opening);

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

        release_domain::<Bn254ScalarField>();

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
        
        //release_domain::<Bn254ScalarField>().unwrap();

        for i in 0..(2 * N) {
            q_g_inv[i] = twist[i].inv() * q_g_inv[i];
        }

        assert_eq!(q_g_inv, target_copy);
    }
}