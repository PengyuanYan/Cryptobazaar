use icicle_core::{msm, msm::MSMConfig};
use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_core::curve::{Curve, Projective, Affine};
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_runtime::memory::{DeviceVec, HostSlice};

use cryptobazaar::utils::load_backend;
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench msmbenchmark */

fn msm_host_benchmark() -> Vec<Projective::<Bn254CurveCfg>> {
    let cfg = MSMConfig::default();
    let mut projective_output = vec![Projective::<Bn254CurveCfg>::zero(); 1];
    let scalars = ScalarCfg::generate_random(10000);
    let bases = Bn254CurveCfg::generate_random_affine_points(10000);
    
    unsafe{
        msm::msm(
            HostSlice::from_slice(&scalars),
            HostSlice::from_slice(&bases),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();
    }

    projective_output
}

fn msm_device_benchmark() -> Vec<Projective::<Bn254CurveCfg>> {
    let cfg = MSMConfig::default();
    let n = 1024usize;
    let mut output_cpu = vec![Projective::<Bn254CurveCfg>::zero(); 1];
    let scalars = ScalarCfg::generate_random(n);
    let bases = Bn254CurveCfg::generate_random_affine_points(1);

    let mut input_scalars_gpu =
        DeviceVec::<Bn254ScalarField>::device_malloc(n).expect("Failed to allocate device memory for input");
    let mut input_bases_gpu =
        DeviceVec::<Affine::<Bn254CurveCfg>>::device_malloc(n).expect("Failed to allocate device memory for input");
    let mut output_gpu =
        DeviceVec::<Projective::<Bn254CurveCfg>>::device_malloc(1).expect("Failed to allocate device memory for output");
    
    input_scalars_gpu
        .copy_from_host(HostSlice::from_slice(&scalars))
        .expect("Failed to copy data to GPU");
    input_bases_gpu
        .copy_from_host(HostSlice::from_slice(&bases))
        .expect("Failed to copy data to GPU");

    unsafe{
        msm::msm(
            &input_scalars_gpu[..],
            &input_bases_gpu[..],
            &cfg,
            &mut output_gpu[..],
        )
        .unwrap();
    }

    output_gpu
        .copy_to_host(HostSlice::from_mut_slice(&mut output_cpu))
        .expect("Failed to copy data back to CPU");

    output_cpu
}

fn criterion_benchmark(c: &mut Criterion) {
    load_backend();

    let id = format!("msmbenchmark");
    c.bench_function(&id, |b| b.iter(|| msm_host_benchmark()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);