use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve, Projective};
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain};
use icicle_runtime::memory::HostSlice;
use icicle_bn254::curve::ScalarCfg;    
use icicle_core::traits::GenerateRandom;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;
use std::ops::Mul;

use cryptobazaar::kzg::{Kzg, PK, VK};
use cryptobazaar::{
    utils::srs::unsafe_setup_from_tau,
    utils::{load_backend, my_ntt, get_device_is_cpu_or_gpu},
};

use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench nonzero */

const N: usize = 128;

fn criterion_benchmark(criterion: &mut Criterion) {
    load_backend();

    let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
    initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
    
    let tau = Bn254ScalarField::from_u32(17u32);
    let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
    let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

    let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
    let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());

    let x = ScalarCfg::generate_random(N);
    let cfg = NTTConfig::<Bn254ScalarField>::default();
    let mut x_coeffs = vec![Bn254ScalarField::zero(); N];
    ntt(
        HostSlice::from_slice(&x),
        NTTDir::kInverse,
        &cfg,
        HostSlice::from_mut_slice(&mut x_coeffs),
    )
    .unwrap();
    let x_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&x_coeffs), x_coeffs.len());

    let id = format!("proof {}", N);
    criterion.bench_function(&id, |b| {
        b.iter(|| {
            /*
                commit to r, 
                commit to r_inv, 
                r_ifft 
                r_coset_fft 
                r_inv_ifft 
                r_inv_coset_fft 
                q_ifft 
                commit to q 
                commit to quotient for kzg opening

                so we can just bench 5 ffts and 4 kzg commits to get realistic bench
             */
            let cpu_or_gpu = get_device_is_cpu_or_gpu();
            let x = ScalarCfg::generate_random(N);
            let cfg = NTTConfig::<Bn254ScalarField>::default();
            let mut x_coeffs = my_ntt::<Bn254CurveCfg>(&x, <Bn254CurveCfg as Curve>::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
            let x_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&x_coeffs), x_coeffs.len());
            let _ = Kzg::commit(&pk, &x_poly);

            let x = ScalarCfg::generate_random(N);
            let mut x_coeffs = my_ntt::<Bn254CurveCfg>(&x, <Bn254CurveCfg as Curve>::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
            let x_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&x_coeffs), x_coeffs.len());
            let _ = Kzg::commit(&pk, &x_poly);

            let x = ScalarCfg::generate_random(N);
            let mut x_coeffs = my_ntt::<Bn254CurveCfg>(&x, <Bn254CurveCfg as Curve>::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
            let x_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&x_coeffs), x_coeffs.len());
            let _ = Kzg::commit(&pk, &x_poly);

            let x = ScalarCfg::generate_random(N);
            let mut x_coeffs = my_ntt::<Bn254CurveCfg>(&x, <Bn254CurveCfg as Curve>::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
            let x_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&x_coeffs), x_coeffs.len());
            let _ = Kzg::commit(&pk, &x_poly);

            let x = ScalarCfg::generate_random(N);
            let mut x_coeffs = my_ntt::<Bn254CurveCfg>(&x, <Bn254CurveCfg as Curve>::ScalarField::one(), NTTDir::kInverse, cpu_or_gpu);
            let x_poly = Bn254Poly::from_coeffs(HostSlice::from_slice(&x_coeffs), x_coeffs.len());
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);