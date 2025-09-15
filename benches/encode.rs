use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve, Projective};
use icicle_bn254::curve::ScalarCfg;    
use icicle_core::traits::GenerateRandom;
use std::ops::Mul;
use icicle_core::{msm, msm::MSMConfig};

use cryptobazaar::utils::load_backend;
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench encode */

fn encode(x: &[Bn254ScalarField], g: &Projective::<Bn254CurveCfg>) -> Vec<Projective::<Bn254CurveCfg>> {
    let mut out: Vec<_> = Vec::with_capacity(x.len());
    for xi in x {
        out.push(g.mul(*xi));
    }
    out
}

fn criterion_benchmark(c: &mut Criterion) {
    load_backend();

    let g = Bn254CurveCfg::get_generator();

    let n_logs = [7, 10, 13];
    for n_logi in n_logs {
        let n = 1 << n_logi;
        let x: Vec<Bn254ScalarField> = ScalarCfg::generate_random(n);

        let id = format!("encode {}", n);
        c.bench_function(&id, |b| b.iter(|| encode(&x, &g)));
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);