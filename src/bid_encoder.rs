use crate::gates::structs::Witness;
use rand::{RngCore, SeedableRng};
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, ntt_inplace, NTT, release_domain};
use icicle_core::curve::{Curve, Affine};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;
use std::marker::PhantomData;

pub struct BidEncoder<const P: usize, const N: usize, C: Curve> {
    pub(crate) bid: [C::ScalarField; N],
    pub(crate) f: [C::ScalarField; N],
    pub(crate) r: [C::ScalarField; N],
}

impl<const P: usize, const N: usize, C: Curve> BidEncoder<P, N, C> {
    pub fn encode<R: RngCore + SeedableRng>(bid: usize, seed: R::Seed) -> Self {
        assert!(bid <= P);
        let mut rng = R::from_seed(seed);
        let length = C::ScalarField::zero().to_bytes_le().len() / 8;  // potential vulnerability

        let mut bid_encoding = vec![C::ScalarField::one(); bid];
        let mut zeroes = vec![C::ScalarField::zero(); P + 1 - bid];
        let mut blinders = Vec::with_capacity(N - P - 1);

        for _ in 0..(N - P - 1) {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            blinders.push(C::ScalarField::from_bytes_le(&bytes));
        }

        bid_encoding.append(&mut zeroes);
        bid_encoding.append(&mut blinders);

        let mut f = Vec::with_capacity(N);
        for _ in 0..N {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            f.push(C::ScalarField::from_bytes_le(&bytes));
        }

        let mut r = Vec::with_capacity(N);
        for _ in 0..N {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            r.push(C::ScalarField::from_bytes_le(&bytes));
        }

        Self {
            bid: bid_encoding.try_into().unwrap(),
            f: f.try_into().unwrap(),
            r: r.try_into().unwrap(),
        }
    }

    pub fn to_gate_witness<R: RngCore + SeedableRng, U>(
        &self,
        seed: R::Seed,
    ) -> Witness<C, U> 
    where
        C: Curve,
        U: UnivariatePolynomial<Field = C::ScalarField>,
        <<C as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C as Curve>::ScalarField>,
        <<C as Curve>::ScalarField as FieldImpl>::Config: NTT<<C as Curve>::ScalarField, <C as Curve>::ScalarField>,
        <C as Curve>::ScalarField: Arithmetic,
    {
        let mut rng = R::from_seed(seed);
        //let domain = get_root_of_unity::<C::ScalarField>((N).try_into().unwrap());
        
        let cfg = NTTConfig::<C::ScalarField>::default();
        
        let mut bid_coeffs = vec![C::ScalarField::zero(); N];
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        ntt(
            HostSlice::from_slice(&self.bid),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut bid_coeffs),
        )
        .unwrap();
        let bid = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&bid_coeffs), bid_coeffs.len());

        let mut f_coeffs = vec![C::ScalarField::zero(); N];
        //domain
        ntt(
            HostSlice::from_slice(&self.f),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut f_coeffs),
        )
        .unwrap();
        let f = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&f_coeffs), f_coeffs.len());


        let mut r_coeffs = vec![C::ScalarField::zero(); N];
        //domain
        ntt(
            HostSlice::from_slice(&self.r),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut r_coeffs),
        )
        .unwrap();
        let r = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&r_coeffs), r_coeffs.len());

        let mut r_inv_evals = self.r[0..P].to_vec();
        for i in 0..P {
            r_inv_evals[i] = r_inv_evals[i].inv();
        }

        let mut r_inv_blinders = Self::sample_blinders(&mut rng, N - P);
        r_inv_evals.append(&mut r_inv_blinders);
        //domain
        ntt_inplace(HostSlice::from_mut_slice(&mut r_inv_evals), NTTDir::kInverse, &cfg).unwrap();
        let r_inv = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&r_inv_evals), r_inv_evals.len());
        
        let mut diff_evals = vec![C::ScalarField::zero(); N];
        let mut g_evals = vec![C::ScalarField::zero(); N];

        for i in 0..P {
            diff_evals[i] = self.bid[i] - self.bid[i + 1];
            g_evals[i] = self.f[i] + self.bid[i] * self.r[i];
        }
        
        let mut diff_blinders = Self::sample_blinders(&mut rng, N - P);
        diff_evals.append(&mut diff_blinders);
        let mut g_blinders = Self::sample_blinders(&mut rng, N - P);
        g_evals.append(&mut g_blinders);

        let mut diff_ntt_evals = vec![C::ScalarField::zero(); N];
        //domain
        ntt(
            HostSlice::from_slice(&diff_evals[..N]),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut diff_ntt_evals),
        )
        .unwrap();

        let mut g_ntt_evals = vec![C::ScalarField::zero(); N];
        //domain
        ntt(
            HostSlice::from_slice(&g_evals[..N]),
            NTTDir::kInverse,
            &cfg,
            HostSlice::from_mut_slice(&mut g_ntt_evals),
        )
        .unwrap();

        //release_domain::<C::ScalarField>().unwrap();

        let diff = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&diff_ntt_evals), diff_ntt_evals.len());
        let g = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&g_ntt_evals), g_ntt_evals.len());

        Witness {
            bid,
            f,
            r,
            r_inv,
            diff,
            g,
            e: PhantomData,
        }
    }

    pub fn to_first_av_round(&self) -> Vec<Affine::<C>> {
        let generator = C::get_generator();
        let mut result = Vec::with_capacity(self.f.len());

        // only sequential
        for fi in &self.f {
            let x = C::mul_scalar(generator, *fi);
            let mut affine_x = Affine::<C>::zero();
            C::to_affine(&x, &mut affine_x);
            result.push(affine_x);
        }
        
        result
    }

    pub fn to_second_av_round(&self, basis: &[Affine::<C>]) -> Vec<Affine::<C>> 
    where
        <C as Curve>::ScalarField: Arithmetic,
    {
        let mut result = Vec::with_capacity(self.f.len());

        for i in 0..self.f.len() {
            let fi = self.f[i];
            let bi = self.bid[i];
            let ri = self.r[i];
            let gi = basis[i];

            let term = fi + bi * ri;
            let x = C::mul_scalar(gi.to_projective(), term);
            let mut affine_x = Affine::<C>::zero();
            C::to_affine(&x, &mut affine_x);
            result.push(affine_x);
        }

        result
    }

    fn sample_blinders<R: RngCore>(rng: &mut R, n: usize) -> Vec<C::ScalarField> {
        let mut v = Vec::with_capacity(n);
        let length = C::ScalarField::zero().to_bytes_le().len() / 8; // potential vulnerability

        for _ in 0..n {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            v.push(C::ScalarField::from_bytes_le(&bytes));
        }

        v
    }
}

#[cfg(test)]
mod encoder_tests {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;

    use rand_chacha::ChaCha20Rng;

    use super::BidEncoder;

    const P: usize = 5;
    const N: usize = 8;

    #[test]
    fn test_encode() {
        let bid = 4usize;
        let enc = [Bn254ScalarField::one(), Bn254ScalarField::one(), Bn254ScalarField::one(), Bn254ScalarField::one(), Bn254ScalarField::zero(), Bn254ScalarField::zero()];

        let seed: [u8; 32] = [
            1, 0, 52, 0, 0, 0, 0, 0, 1, 0, 10, 0, 22, 32, 0, 0, 2, 0, 55, 49, 0, 11, 0, 0, 3, 0, 0,
            0, 0, 0, 2, 92,
        ];

        let encoder = BidEncoder::<P, N, Bn254CurveCfg>::encode::<ChaCha20Rng>(bid, seed);
        assert_eq!(enc, encoder.bid[0..(P + 1)]);
    }
}