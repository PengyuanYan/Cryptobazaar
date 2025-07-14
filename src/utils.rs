use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::curve::Curve;
use icicle_core::ntt::{NTTDomain, get_root_of_unity};
use icicle_core::traits::Arithmetic;
use std::ops::Mul;

pub mod srs;

pub fn is_pow_2(x: usize) -> bool {
    (x & (x - 1)) == 0
}

pub fn evaluate_vanishing_over_extended_coset<C: Curve>(n: usize, k: usize) -> Vec<C::ScalarField> 
where 
    <C as Curve>::ScalarField: Arithmetic,
    <<C as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C as Curve>::ScalarField>,
{
    assert!(is_pow_2(n));
    assert!(is_pow_2(k));

    let p = C::ScalarField::from_u32(5u32);
    let domain_kn = get_root_of_unity::<C::ScalarField>((k * n).try_into().unwrap());

    let coset_generator_pow_n = p.pow(n.try_into().unwrap());
    let wi = domain_kn;

    let mut modulus_zh_coset_evals = Vec::with_capacity(k);

    for i in 0usize..k {
        let zhi = coset_generator_pow_n * wi.pow((i * n).try_into().unwrap()) - C::ScalarField::one();
        modulus_zh_coset_evals.push(zhi);
    }

    modulus_zh_coset_evals
}



pub fn get_coeffs_of_poly<P>(poly: &P) -> Vec<P::Field>
where
    P: UnivariatePolynomial,
    P::Field: FieldImpl,
{
    let n: usize = poly.get_nof_coeffs().try_into().expect("Host's archtecture is smaller than 64-bit");

    let mut coeffs = vec![P::Field::zero(); n];

    poly.copy_coeffs(0, HostSlice::from_mut_slice(&mut coeffs));

    coeffs
}