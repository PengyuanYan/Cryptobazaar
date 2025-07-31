use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve, Affine, Projective};
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;
use icicle_bn254::curve::ScalarCfg;    
use icicle_core::traits::GenerateRandom;
use std::ops::Mul;


use criterion::{criterion_group, criterion_main, Criterion};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* RUN WITH: cargo bench --bench encode */

fn encode(x: &[F], g: &G1) -> Vec<G1> {
    let mut out: Vec<_> = Vec::with_capacity(x.len());
    for xi in &x {
        out.push(g.mul(xi));
    }
    out
}

fn criterion_benchmark(c: &mut Criterion) {
    let g = Bn254G2CurveCfg::get_generator();

    let n_logs = [7, 10, 13];
    for n_logi in n_logs {
        let n = 1 << n_logi;
        let x: Vec<F> = ScalarCfg::generate_random(n);

        let id = format!("encode {}", n);
        c.bench_function(&id, |b| b.iter(|| encode(&x, &g)));
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);