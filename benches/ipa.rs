use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;
use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;
use icicle_bn254::curve::ScalarCfg;
use icicle_core::traits::GenerateRandom;
use icicle_core::{msm, msm::MSMConfig};
use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, release_domain, NTT};
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::pairing::Pairing;


use cryptobazaar::{
    ipa::{
        structs::{Instance, Witness},
        InnerProduct,
    },
    kzg::PK,
    utils::srs::unsafe_setup_from_tau,
    utils::{evaluate_all_lagrange_coefficients, load_backend, get_device_is_cpu_or_gpu},
};
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench ipa */

const N: usize = 8192;
const LOG_N: usize = 13;

fn prove<const N: usize, C1, C2, F, U>(
    instance: &Instance<N, C1>,
    witness: &Witness<N, C1>,
    pk: &PK<C1, C2, F>,
    cpu_or_gpu: usize
) 
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
    <C1 as Curve>::ScalarField: Arithmetic,
    <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    <U as UnivariatePolynomial>::FieldConfig: GenerateRandom<<C1 as Curve>::ScalarField>,
    <<C1 as Curve>::ScalarField as FieldImpl>::Config: GenerateRandom<<C1 as Curve>::ScalarField>
{
    InnerProduct::<N, LOG_N, C1, C2, F, U>::prove(instance, witness, pk, cpu_or_gpu);
}

fn criterion_benchmark(criterion: &mut Criterion) {
    load_backend();
    let cpu_or_gpu = get_device_is_cpu_or_gpu();

    let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
    initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

    let domain = get_root_of_unity::<Bn254ScalarField>((N).try_into().unwrap());

    let tau = Bn254ScalarField::from_u32(100u32);
    let lb_at_tau = evaluate_all_lagrange_coefficients::<Bn254CurveCfg>(domain, tau, N);

    let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);

    let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };

    let generator = Bn254CurveCfg::get_generator();
    let a: Vec<Bn254ScalarField> = ScalarCfg::generate_random(N);

    let mut lagrange_basis: Vec<Affine<Bn254CurveCfg>> = Vec::with_capacity(lb_at_tau.len());
    for li in lb_at_tau {
        let point = Bn254CurveCfg::mul_scalar(generator, li);
        lagrange_basis.push(point.into());
    }

    let b = Bn254CurveCfg::generate_random_affine_points(N);
    let h_base_projective = Bn254CurveCfg::mul_scalar(generator, ScalarCfg::generate_random(1)[0]);
    let mut h_base_affine = Affine::<Bn254CurveCfg>::zero();
    Bn254CurveCfg::to_affine(&h_base_projective, &mut h_base_affine);

    let cfg = MSMConfig::default();
    let mut ac_projective = vec![Projective::<Bn254CurveCfg>::zero(); 1];    
    msm::msm(
        HostSlice::from_slice(&a),
        HostSlice::from_slice(&lagrange_basis),
        &cfg,
        HostSlice::from_mut_slice(&mut ac_projective),
    )
    .unwrap();
    
    let mut ac_affine = Affine::<Bn254CurveCfg>::zero();
    Bn254CurveCfg::to_affine(&ac_projective[0], &mut ac_affine);

    let mut c: Vec<Affine::<Bn254CurveCfg>> = Vec::with_capacity(b.len());
    for i in 0..b.len() {
        let bi = b[i];
        let ai = a[i];

        let c_projective = Bn254CurveCfg::mul_scalar(bi.to_projective(), ai);
        let mut c_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&c_projective, &mut c_affine);
        c.push(c_affine);
    }

    let instance = Instance::<N, Bn254CurveCfg> {
        ac: ac_affine,
        b: b.try_into().unwrap(),
        h_base: h_base_affine,
        c: c.try_into().unwrap(),
    };

    let witness = Witness::<N, Bn254CurveCfg> {
        a: a.try_into().unwrap(),
    };

    let id = format!("proof {}", N);
    criterion.bench_function(&id, |b| {
        b.iter(|| prove::<N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>(&instance, &witness, &pk, cpu_or_gpu))
    });

    release_domain::<Bn254ScalarField>().unwrap();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);