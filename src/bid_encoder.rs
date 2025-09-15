//! Bid_encoder module is the front door module of other module.
//! It encode the bid values to bid vector and sample the random vector x and r.
//!
//! This module provides a encoding interface with:
//! - `encode` — deterministically encode the bid values to bid vector and sample the random vector x and r
//! - `to_gate_witness` — further encode the bid vector, the random vector as the gate witness
//! - `to_first_av_round` — finish the duty of bidders for the first round anonymous veto
//! - `to_second_av_round` — finish the duty of bidders for the second round anonymous veto
//! - `sample_blinders` — proved blinding vectors
use crate::gates::structs::Witness;
use rand::{RngCore, SeedableRng};
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, ntt_inplace, NTT, release_domain};
use icicle_core::curve::{Curve, Affine};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;
use std::marker::PhantomData;
use crate::utils::{my_ntt, get_device_is_cpu_or_gpu};

pub struct BidEncoder<const P: usize, const N: usize, C: Curve> {
    pub(crate) bid: [C::ScalarField; N],
    pub(crate) random_x: [C::ScalarField; N],
    pub(crate) random_r: [C::ScalarField; N],
}

impl<const P: usize, const N: usize, C: Curve> BidEncoder<P, N, C> {
    /// Encode bid deterinistically according to a give seed.
    pub fn encode<R: RngCore + SeedableRng>(bid: usize, seed: R::Seed) -> Self {
        assert!(bid <= P);
        let mut rng = R::from_seed(seed);
        let length = C::ScalarField::zero().to_bytes_le().len() / 8;  // potential vulnerability
        
        // construct an unary bid
        let mut bid_encoding = vec![C::ScalarField::one(); bid];
        let mut zeroes = vec![C::ScalarField::zero(); P + 1 - bid];
        let mut blinders = Vec::with_capacity(N - P - 1);

        for _ in 0..(N - P - 1) {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            blinders.push(C::ScalarField::from_bytes_le(&bytes));
        }
        
        bid_encoding.append(&mut zeroes);
        // blind the bid vecotr
        bid_encoding.append(&mut blinders);
        
        // generate random vector x
        let mut random_x = Vec::with_capacity(N);
        for _ in 0..N {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            random_x.push(C::ScalarField::from_bytes_le(&bytes));
        }
        
        // generate random vector r
        let mut random_r = Vec::with_capacity(N);
        for _ in 0..N {
            let mut bytes = vec![0u8; length * 4];
            rng.fill_bytes(&mut bytes);
            random_r.push(C::ScalarField::from_bytes_le(&bytes));
        }

        Self {
            bid: bid_encoding.try_into().unwrap(),
            random_x: random_x.try_into().unwrap(),
            random_r: random_r.try_into().unwrap(),
        }
    }
    
    /// Further deterministically encode the bid and random vectors.
    /// Make them is ready for the prove function of gate.
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
        let cpu_or_gpu = get_device_is_cpu_or_gpu();

        let mut rng = R::from_seed(seed);
        
        //domain
        let bid_coeffs = my_ntt::<C>(&self.bid, C::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let bid = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&bid_coeffs), bid_coeffs.len());
        
        //domain
        let random_x_coeffs = my_ntt::<C>(&self.random_x, C::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let random_x = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&random_x_coeffs), random_x_coeffs.len());
        
        //domain
        let random_r_coeffs = my_ntt::<C>(&self.random_r, C::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let random_r = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&random_r_coeffs), random_r_coeffs.len());

        let mut random_r_inv_evals = self.random_r[0..P].to_vec();
        for i in 0..P {
            random_r_inv_evals[i] = random_r_inv_evals[i].inv();
        }

        let mut random_r_inv_blinders = Self::sample_blinders(&mut rng, N - P);
        random_r_inv_evals.append(&mut random_r_inv_blinders);
        
        //domain
        let random_r_inv_evals = my_ntt::<C>(&random_r_inv_evals, C::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let random_r_inv = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&random_r_inv_evals), random_r_inv_evals.len());
        
        let mut diff_f_evals = vec![C::ScalarField::zero(); N];
        let mut hidden_bid_evals = vec![C::ScalarField::zero(); N];

        // diff_f_evals[i] = bid[i] - bid[i+1]
        diff_f_evals
            .iter_mut()
            .zip(self.bid.windows(2).take(N))
            .for_each(|(d, w)| *d = w[0] - w[1]);

        // hidden_bid_evals[i] = random_x[i] + bid[i] * random_r[i]
        hidden_bid_evals
            .iter_mut()
            .zip(self.random_x.iter()
                 .zip(self.bid.iter().zip(self.random_r.iter()))
                 .take(N))
            .for_each(|(out, (x, (b, r)))| *out = *x + (*b) * (*r));
        
        let mut diff_f_blinders = Self::sample_blinders(&mut rng, N - P);
        diff_f_evals.append(&mut diff_f_blinders);
        let mut g_blinders = Self::sample_blinders(&mut rng, N - P);
        hidden_bid_evals.append(&mut g_blinders);
        
        //domain
        let diff_f_ntt_evals = my_ntt::<C>(&diff_f_evals[..N], C::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let diff_f = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&diff_f_ntt_evals), diff_f_ntt_evals.len());
        
        //domain
        let hidden_bid_ntt_evals = my_ntt::<C>(&hidden_bid_evals[..N], C::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
        let hidden_bid = UnivariatePolynomial::from_coeffs(HostSlice::from_slice(&hidden_bid_ntt_evals), hidden_bid_ntt_evals.len());

        Witness {
            bid,
            random_x,
            random_r,
            random_r_inv,
            diff_f,
            hidden_bid,
            e: PhantomData,
        }
    }
    
    /// Compute the X for the first round anonymous veto.
    pub fn to_first_av_round(&self) -> Vec<Affine::<C>> {
        let generator = C::get_generator();
        let mut result = Vec::with_capacity(self.random_x.len());

        // only sequential
        for fi in &self.random_x {
            let x = C::mul_scalar(generator, *fi);
            let mut affine_x = Affine::<C>::zero();
            C::to_affine(&x, &mut affine_x);
            result.push(affine_x);
        }
        
        result
    }
    
    /// Compute the Z for the second round anonymous veto.
    pub fn to_second_av_round(&self, basis: &[Affine::<C>]) -> Vec<Affine::<C>> 
    where
        <C as Curve>::ScalarField: Arithmetic,
    {
        let mut result = Vec::with_capacity(self.random_x.len());

        for i in 0..self.random_x.len() {
            let xi = self.random_x[i];
            let bi = self.bid[i];
            let ri = self.random_r[i];
            let yi = basis[i];
            
            // hidden_bid
            let term = xi + bi * ri;
            // Z = hidden_bid * Y
            let x = C::mul_scalar(yi.to_projective(), term);
            let mut affine_x = Affine::<C>::zero();
            C::to_affine(&x, &mut affine_x);
            result.push(affine_x);
        }

        result
    }
    
    /// Deterministically ample a randome vector according to the input seed.
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
    use icicle_core::curve::{Curve, Affine};
    use icicle_core::traits::FieldImpl;
    use rand_chacha::ChaCha20Rng;

    use super::BidEncoder;

    const P: usize = 8;
    const N: usize = 12;

    fn seed() -> [u8; 32] {
        [
            1, 0, 77, 0, 0, 1, 0, 0, 0, 0, 10, 0, 12, 32, 0, 0, 2, 0, 55, 49, 0, 11, 0, 0, 1, 9, 1,
            0, 0, 0, 2, 6,
        ]
    }

    #[test]
    fn test_encode_layout_prefix() {
        let bid = 5usize;
        let encoder = BidEncoder::<P, N, Bn254CurveCfg>::encode::<ChaCha20Rng>(bid, seed());

        let expected_prefix = [
            Bn254ScalarField::one(),
            Bn254ScalarField::one(),
            Bn254ScalarField::one(),
            Bn254ScalarField::one(),
            Bn254ScalarField::one(),

            Bn254ScalarField::zero(),
            Bn254ScalarField::zero(),
            Bn254ScalarField::zero(),
            Bn254ScalarField::zero(),
        ];

        assert_eq!(&expected_prefix, &encoder.bid[0..(P + 1)]);
        assert_eq!(encoder.bid.len(), N);
        assert_eq!(encoder.random_x.len(), N);
        assert_eq!(encoder.random_r.len(), N);
    }
}
