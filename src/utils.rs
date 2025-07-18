use icicle_runtime::memory::HostSlice;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::curve::Curve;
use icicle_core::ntt::{NTTDomain, get_root_of_unity};
use icicle_core::traits::Arithmetic;
use std::ops::{Mul, Sub, Add};
use std::slice::IterMut;

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


/// Evaluate all the lagrange polynomials defined by this domain at the
    /// point `tau`. This is computed in time O(|domain|).
    /// Then given the evaluations of a degree d polynomial P over this domain,
    /// where d < |domain|, `P(tau)` can be computed as
    /// `P(tau) = sum_{i in [|Domain|]} L_{i, Domain}(tau) * P(g^i)`.
    /// `L_{i, Domain}` is the value of the i-th lagrange coefficient
    /// in the returned vector.
pub fn evaluate_all_lagrange_coefficients<C: Curve>(domain: C::ScalarField, tau: C::ScalarField, n: usize) -> Vec<C::ScalarField>
where
    <C as Curve>::ScalarField: Arithmetic,
    C::ScalarField: Mul<Output = C::ScalarField>,
    C::ScalarField: Sub<Output = C::ScalarField>,
    C::ScalarField: Add<Output = C::ScalarField>
{
    // Evaluate all Lagrange polynomials at tau to get the lagrange coefficients.
        // Define the following as
        // - H: The coset we are in, with generator g and offset h
        // - m: The size of the coset H
        // - Z_H: The vanishing polynomial for H. Z_H(x) = prod_{i in m} (x - hg^i) = x^m - h^m
        // - v_i: A sequence of values, where v_0 = 1/(m * h^(m-1)), and v_{i + 1} = g * v_i
        //
        // We then compute L_{i,H}(tau) as `L_{i,H}(tau) = Z_H(tau) * v_i / (tau - h * g^i)`
        //
        // However, if tau in H, both the numerator and denominator equal 0
        // when i corresponds to the value tau equals, and the coefficient is 0
        // everywhere else. We handle this case separately, and we can easily
        // detect by checking if the vanishing poly is 0.
    let size = n;
    let z_h_at_tau = tau.pow(n) - C::ScalarField::one();
    let offset = C::ScalarField::one();
    let group_gen = domain;

    assert_eq!(domain.pow(n), C::ScalarField::one());

    if z_h_at_tau == C::ScalarField::zero() {
        // In this case, we know that tau = hg^i, for some value i.
            // Then i-th lagrange coefficient in this case is then simply 1,
            // and all other lagrange coefficients are 0.
            // Thus we find i by brute force.
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
        // In this case we have to compute `Z_H(tau) * v_i / (tau - h g^i)`
            // for i in 0..size
            // We actually compute this by computing (Z_H(tau) * v_i)^{-1} * (tau - h g^i)
            // and then batch inverting to get the correct lagrange coefficients.
            // We let `l_i = (Z_H(tau) * v_i)^-1` and `r_i = tau - h g^i`
            // Notice that since Z_H(tau) is i-independent,
            // and v_i = g * v_{i-1}, it follows that
            // l_i = g^-1 * l_{i-1}
            // TODO: consider caching the computation of l_i to save N multiplications
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