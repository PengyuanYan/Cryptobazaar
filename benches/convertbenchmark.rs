use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;
use icicle_runtime::memory::DeviceSlice;
use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_core::curve::{Curve, Projective, Affine};
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;
use std::ops::Mul;

use cryptobazaar::utils::load_backend;
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench convertbenchmark */

fn criterion_benchmark(c: &mut Criterion) {
    let bases_1 = Bn254CurveCfg::generate_random_affine_points(8192 + 8192);
    let bases_2 = Bn254CurveCfg::generate_random_affine_points(8192 + 8192);
    let scalars = ScalarCfg::generate_random(8192 + 8192);
    let mut bases = Vec::<Projective::<Bn254CurveCfg>>::with_capacity(8192 + 8192);
    
    for i in 0..(8192 + 8192) {
        bases.push(bases_1[i].to_projective() + bases_2[i].to_projective());
    }

    for i in 0..(8192 + 8192) {
        bases[i] = bases[i] + bases_1[i].to_projective().mul(scalars[i]);
    }

    let id = format!("msmbenchmark");
    c.bench_function(&id, |b| {
        b.iter(|| {
            let mut affine = Vec::<Affine::<Bn254CurveCfg>>::with_capacity(8192 + 8192);
            for i in 0..(8192) {
                bases[i] = bases[i] + bases_1[i].to_projective().mul(scalars[i]);
                let base: Projective::<Bn254CurveCfg> = bases[i].into();
                let mut affine_point = Affine::<Bn254CurveCfg>::zero();
                Bn254CurveCfg::to_affine(&base ,&mut affine_point);
                affine.push(affine_point);
            }
    })});
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
