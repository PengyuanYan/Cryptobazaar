use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_core::curve::{Curve, Projective, Affine};
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_runtime::memory::{DeviceVec, HostSlice};

use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, release_domain};
use icicle_core::traits::FieldImpl;

use cryptobazaar::utils::load_backend;
use criterion::{criterion_group, criterion_main, Criterion};

const N: usize = 8192;

/* RUN WITH: cargo bench --bench nttbenchmark */

fn ntt_host_benchmark() -> Vec<Bn254ScalarField> {
    let scalars = ScalarCfg::generate_random(N);

    let cfg = NTTConfig::<Bn254ScalarField>::default();
    let mut coeffs = vec![Bn254ScalarField::zero(); N];
    ntt(
        HostSlice::from_slice(&scalars),
        NTTDir::kForward,
        &cfg,
        HostSlice::from_mut_slice(&mut coeffs),
    )
    .unwrap();

    coeffs
}

fn ntt_device_benchmark() -> Vec<Bn254ScalarField> {
    let scalars = ScalarCfg::generate_random(N);

    let cfg = NTTConfig::<Bn254ScalarField>::default();

    let mut input_gpu =
        DeviceVec::<Bn254ScalarField>::device_malloc(N).expect("Failed to allocate device memory for ntt input");
    let mut output_gpu =
        DeviceVec::<Bn254ScalarField>::device_malloc(N).expect("Failed to allocate device memory for ntt output");
    input_gpu
        .copy_from_host(HostSlice::from_slice(&scalars))
        .expect("Failed to copy ntt input data to GPU");
    
    let mut coeffs = vec![Bn254ScalarField::zero(); N];
    ntt(
        &input_gpu[..],
        NTTDir::kForward,
        &cfg,
        &mut output_gpu[..],
    )
    .unwrap();

    output_gpu
        .copy_to_host(HostSlice::from_mut_slice(&mut coeffs))
        .expect("Failed to copy ntt result data back to CPU");

    coeffs
}

fn criterion_benchmark(c: &mut Criterion) {
    load_backend();
    
    let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
    initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
    
    let id = format!("nttbenchmark");
    c.bench_function(&id, |b| b.iter(|| ntt_host_benchmark()));

    release_domain::<Bn254ScalarField>().unwrap();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);