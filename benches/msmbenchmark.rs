use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;
use icicle_runtime::memory::DeviceSlice;
use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_core::curve::{Curve, Projective};
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;

use cryptobazaar::utils::load_backend;
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench msmbenchmark */

fn msm_host_benchmark() -> Vec<Projective::<Bn254CurveCfg>> {
    let cfg = MSMConfig::default();
    let mut projective_output = vec![Projective::<Bn254CurveCfg>::zero(); 1];
    let scalars = ScalarCfg::generate_random(1000);
    let bases = Bn254CurveCfg::generate_random_affine_points(1000);
    
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
    let mut projective_output = vec![Projective::<Bn254CurveCfg>::zero(); 1];
    let scalars = ScalarCfg::generate_random(1000);
    let bases = Bn254CurveCfg::generate_random_affine_points(1000);
    
    unsafe{
        msm::msm(
            DeviceSlice::from_slice(&scalars),
            DeviceSlice::from_slice(&bases),
            &cfg,
            DeviceSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();
    }

    projective_output
}

fn criterion_benchmark(c: &mut Criterion) {
    load_backend();

    let id = format!("msmbenchmark");
    c.bench_function(&id, |b| b.iter(|| msm_host_benchmark()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);