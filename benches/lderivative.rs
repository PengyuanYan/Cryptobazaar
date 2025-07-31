use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
use icicle_bn254::curve::ScalarField as Bn254ScalarField;
use icicle_core::curve::Curve;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, release_domain, NTT};
use icicle_runtime::memory::HostSlice;
use icicle_core::traits::Arithmetic;
use icicle_core::pairing::Pairing;
// use icicle_core::{msm, msm::MSMConfig};

use cryptobazaar::kzg::{Kzg, PK, VK};
use cryptobazaar::{
    zk_log_derivative::{
        structs::{Instance, Witness, ProverIndex, VerifierIndex},
        Argument,
    },
    utils::srs::unsafe_setup_from_tau,
    utils::load_backend,
};
use criterion::{criterion_group, criterion_main, Criterion};

/* RUN WITH: cargo bench --bench lderivative */

const N: usize = 8192;
const B: usize = 1;

fn prove<const N: usize, C1, C2, F, P>(
    index_p: &ProverIndex::<C1, P>,
    index_v: &VerifierIndex<C1>,
    instance: &Instance<C1>,
    witness: &Witness<C1, P>,
    pk: &PK<C1, C2, F>,
)
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    P: UnivariatePolynomial<Field = C1::ScalarField>,
    <C1 as Curve>::ScalarField: Arithmetic,
    <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>
{
    let _ = Argument::<N, B, C1, C2, F, P>::prove(index_p, index_v, &instance, &witness, &pk);
}

fn criterion_benchmark(criterion: &mut Criterion) {
    load_backend();

    let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
    initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

    let domain = get_root_of_unity::<Bn254ScalarField>(N.try_into().unwrap());

    let tau = Bn254ScalarField::from_u32(17u32);
    let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
    let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

    let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData,};
    let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());

    let index_v = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::index_v(&pk);
    let index_p = Argument::<N, B, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::index_p();

    let mut f_evals = vec![Bn254ScalarField::zero(); N - B];
    f_evals[3] = Bn254ScalarField::one();

    let mut blinders: Vec<_> = (0..B).map(|i| Bn254ScalarField::from_u32((i + 10) as u32)).collect();
    let blinders_cloned = blinders.clone();
    f_evals.append(&mut blinders);
    
    let cfg = NTTConfig::<Bn254ScalarField>::default();
    let mut f_coeffs = vec![Bn254ScalarField::zero(); N];
    ntt(
        HostSlice::from_slice(&f_evals),
        NTTDir::kInverse,
        &cfg,
        HostSlice::from_mut_slice(&mut f_coeffs),
    )
    .unwrap();

    let f = Bn254Poly::from_coeffs(HostSlice::from_slice(&f_coeffs), N);
    let f_cm = Kzg::commit(&pk, &f);

    let instance = Instance::<Bn254CurveCfg> { f_cm };

    let witness = Witness::<Bn254CurveCfg, Bn254Poly> { f, e: PhantomData,};

    let id = format!("proof {}", N);
    criterion.bench_function(&id, |b| {
        b.iter(|| prove::<N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>(&index_p, &index_v, &instance, &witness, &pk))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);