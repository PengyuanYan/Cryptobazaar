use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve, Affine};
use icicle_core::traits::FieldImpl;
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;
use std::ops::Mul;

use cryptobazaar::auctioneer::Auctioneer;
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench auctioneer_r1 */

fn setup_round_1<const N: usize, const B: usize>() -> Auctioneer<N, B, Bn254CurveCfg> {
    let g = Bn254CurveCfg::get_generator();

    let mut a = Auctioneer::<N, B, Bn254CurveCfg>::new();
    let mut secrets = vec![vec![Bn254ScalarField::zero(); N]; B];
    let mut first_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];

    // initialize n msgs fro each party
    for i in 0..B {
        for j in 0..N {
            secrets[i][j] = ScalarCfg::generate_random(1)[0];
        }
    }

    // initialize n msgs fro each party
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

fn setup_round_2<const N: usize, const B: usize>() -> Auctioneer<N, B, Bn254CurveCfg> {
    let g = Bn254CurveCfg::get_generator();

    let mut a = Auctioneer::<N, B, Bn254CurveCfg>::new();
    let mut secrets = vec![vec![Bn254ScalarField::zero(); N]; B];
    let mut first_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];

    // initialize n msgs fro each party
    for i in 0..B {
        for j in 0..N {
            secrets[i][j] = ScalarCfg::generate_random(1)[0];
        }
    }

    // initialize n msgs fro each party
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

    // we get output for each party per round
    // where each row is of len B (output of av for each party)
    let fr_result = a.output_first_round();

    let mut second_msgs = vec![vec![Affine::<Bn254CurveCfg>::zero(); N]; B];
    for i in 0..B {
        for j in 0..N {
            let projective_result = fr_result[j][i].to_projective().mul(secrets[i][j]);
            let mut affine_result = Affine::<Bn254CurveCfg>::zero();
            Bn254CurveCfg::to_affine(&projective_result, &mut affine_result);
            second_msgs[i][j] = affine_result;
        }
    }

    // each party sends it's second round msgs
    for i in 0..B {
        a.register_msgs(&second_msgs[i], i).unwrap();
    }

    a
}

fn bench_second_round<const N: usize, const B: usize>(
    a: Auctioneer<N, B, Bn254CurveCfg>,
) -> Vec<Affine::<Bn254CurveCfg>> {
    let mut a_clone = a.clone();
    a_clone.output_second_round()
}

fn bench_first_round<const N: usize, const B: usize>(
    a: Auctioneer<N, B, Bn254CurveCfg>,
) -> Vec<Vec<Affine::<Bn254CurveCfg>>> {
    let mut a_clone = a.clone();
    a_clone.output_first_round()
}

fn round_1(c: &mut Criterion) {
    const N: usize = 8192;
    const B: usize = 256;

    let a = setup_round_1::<N, B>();
    let id = format!("Round1: range = {}, bidders = {}", N, B);
    c.bench_function(&id, |b| b.iter(|| bench_first_round(a.clone())));
}

fn round_2(c: &mut Criterion) {
    const N: usize = 32;
    const B: usize = 32;

    let a = setup_round_2::<N, B>();
    let id = format!("Round2: range = {}, bidders = {}", N, B);
    c.bench_function(&id, |b| b.iter(|| bench_second_round(a.clone())));
}

fn criterion_benchmark(c: &mut Criterion) {
    //round_1(c);
    round_2(c);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);