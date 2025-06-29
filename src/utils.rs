use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;

pub mod srs;

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