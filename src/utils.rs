use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::curve::Curve;
use icicle_core::ntt::{NTTDomain, get_root_of_unity};
use icicle_core::traits::Arithmetic;

pub mod folding;
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


// this function is reimplemented based on the source code of ark library
pub fn evaluate_all_lagrange_coefficients<C: Curve>(domain: C::ScalarField, tau: C::ScalarField, n: usize) -> Vec<C::ScalarField>
where
    <C as Curve>::ScalarField: Arithmetic,
{
    let size = n;
    let z_h_at_tau = tau.pow(n) - C::ScalarField::one();
    let offset = C::ScalarField::one();
    let group_gen = domain;

    assert_eq!(domain.pow(n), C::ScalarField::one());

    if z_h_at_tau == C::ScalarField::zero() {
        let mut u = vec![C::ScalarField::zero(); size];
        let mut omega_i = offset;
        for u_i in u.iter_mut().take(size) {
            if omega_i == tau {
                *u_i = C::ScalarField::one();
                break;
            }
            omega_i = omega_i * group_gen;
        }
        
        u

    } else {
        let group_gen_inv = domain.inv();

        // v_0_inv = m * h^(m-1)
        let v_0_inv = C::ScalarField::from_u32(n as u32) * offset.pow(size - 1);
        let mut l_i = z_h_at_tau.inv() * v_0_inv;
        let mut negative_cur_elem = C::ScalarField::zero() - offset;
        let mut lagrange_coefficients_inverse = vec![C::ScalarField::zero(); size];
        for coeff in &mut lagrange_coefficients_inverse {
            let r_i = tau + negative_cur_elem;
            *coeff = l_i * r_i;
            // Increment l_i and negative_cur_elem
            l_i = l_i * group_gen_inv;
            negative_cur_elem = negative_cur_elem * group_gen;
        }

        // Invert the lagrange coefficients inverse, to get the actual coefficients,
        // and return these
        for i in 0..lagrange_coefficients_inverse.len() {
            lagrange_coefficients_inverse[i] = lagrange_coefficients_inverse[i].inv();
        }

        lagrange_coefficients_inverse
    }
}