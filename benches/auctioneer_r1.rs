use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve, Affine, Projective};
use icicle_core::traits::FieldImpl;
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;
use std::ops::Mul;

use cryptobazaar::auctioneer::Auctioneer;
use criterion::{criterion_group, criterion_main, Criterion};

use cryptobazaar::utils::load_backend;

/* RUN WITH: cargo bench --bench auctioneer_r1 */

fn setup_round_1<const N: usize, const B: usize>() -> Auctioneer<N, B, Bn254CurveCfg> {
    let g = Bn254CurveCfg::get_generator();

    let mut a = Auctioneer::<N, B, Bn254CurveCfg>::new();
    let mut secrets = vec![vec![Bn254ScalarField::zero(); N]; B];
    let mut first_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];
    
    // initialize n msgs for each party
    for i in 0..B {
        let generated_secrets = ScalarCfg::generate_random(N);
        for j in 0..N {
            secrets[i][j] = generated_secrets[j];
        }
    }

    // initialize n msgs for each party
    for i in 0..B {
        for j in 0..N {
            let projective_result = g.mul(secrets[i][j]);
            let mut affine_result = Affine::<Bn254CurveCfg>::zero();
            Bn254CurveCfg::to_affine(&projective_result, &mut affine_result);
            first_msgs[i][j] = affine_result;
        }
    }

    // each party sends it's first round msgs
    for i in 0..B {
        a.register_msgs(&first_msgs[i], i).unwrap();
    }

    a
}

fn bench_first_round<const N: usize, const B: usize>(
    a: Auctioneer<N, B, Bn254CurveCfg>,
) -> Vec<Vec<Affine::<Bn254CurveCfg>>> {
    let mut a_clone = a.clone();
    a_clone.output_first_round()
}

fn round_1(c: &mut Criterion) {
    const N: usize = 1024;
    const B: usize = 128;

    let a = setup_round_1::<N, B>();
    let id = format!("Round1: range = {}, bidders = {}", N, B);
    c.bench_function(&id, |b| b.iter(|| bench_first_round(a.clone())));
}

fn criterion_benchmark(c: &mut Criterion) {
    load_backend();
    round_1(c);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);