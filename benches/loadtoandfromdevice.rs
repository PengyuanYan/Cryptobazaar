use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve, Projective, Affine};
use icicle_bn254::curve::ScalarCfg;    
use icicle_core::traits::GenerateRandom;
use std::ops::Mul;
use std::time::{Duration, Instant};
use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::{DeviceVec, HostSlice};

use cryptobazaar::utils::{load_backend, get_device_is_cpu_or_gpu};
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench loadtoandfromdevice */

const N: usize = 8192;

fn criterion_benchmark(c: &mut Criterion) {
    load_backend();
    if get_device_is_cpu_or_gpu() == 0 {
        panic!("this test need gpu");
    }

    let scalars: Vec<Bn254ScalarField> = ScalarCfg::generate_random(N);
    let bases: Vec<Affine::<Bn254CurveCfg>> = Bn254CurveCfg::generate_random_affine_points(N);
    
    let id = format!("Size = {}", N);
    c.bench_function(&id, |b| {
        b.iter_custom(|iters| {
            let mut total = Duration::ZERO;
            let start = Instant::now();
            let cfg = MSMConfig::default();
            let mut output_cpu = vec![Projective::<Bn254CurveCfg>::zero(); 1];
            let mut input_scalars_gpu =
                DeviceVec::<<Bn254CurveCfg as Curve>::ScalarField>::device_malloc(N).unwrap();
            let mut input_bases_gpu =
                DeviceVec::<Affine::<Bn254CurveCfg>>::device_malloc(N).unwrap();
            let mut output_gpu =
                DeviceVec::<Projective::<Bn254CurveCfg>>::device_malloc(1).unwrap();
        
            input_scalars_gpu
                .copy_from_host(HostSlice::from_slice(&scalars)).unwrap();
            input_bases_gpu
                .copy_from_host(HostSlice::from_slice(&bases)).unwrap();
            total += start.elapsed();

            msm::msm(
                &input_scalars_gpu[..],
                &input_bases_gpu[..],
                &cfg,
                &mut output_gpu[..],
            )
            .unwrap();
        
            let start = Instant::now();
            output_gpu
                .copy_to_host(HostSlice::from_mut_slice(&mut output_cpu)).unwrap();
            total += start.elapsed();
        
            total
        })
    });

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);