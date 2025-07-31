use icicle_core::curve::{Curve, Affine, Projective};
use icicle_core::traits::FieldImpl;
use icicle_core::ntt::{NTTDomain, get_root_of_unity};
use icicle_core::traits::Arithmetic;
use crate::utils::is_pow_2;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_runtime::memory::HostSlice;

pub fn compute_lagrange_basis_commitments<C>(tau_powers: &[Projective::<C>], ) -> Vec<Affine::<C>>
where
    C: Curve,
    C::ScalarField: FieldImpl,
    <C as Curve>::ScalarField: Arithmetic,
    <<C as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C as Curve>::ScalarField>,
{
    let n = tau_powers.len();
    assert!(is_pow_2(n));

    let omega = get_root_of_unity::<C::ScalarField>((n).try_into().unwrap());
    let psi = omega.inv();
    let n_inv = C::ScalarField::from_u32((n).try_into().unwrap()).inv();
    
    let mut out = Vec::with_capacity(n);
    
    // L_i(τ)G = (1/n) * Σ_j ψ^{i*j} · [τ^j]G
    for i in 0..n {
        // r = ψ^i
        let r = psi.pow((i).try_into().unwrap());
        let mut w = C::ScalarField::one();

        let mut acc = Projective::<C>::zero();
        for p_j in tau_powers {
            acc = acc + C::mul_scalar(*p_j, w);
            w = w * r;
        }
        
        // scale by 1/n
        acc = C::mul_scalar(acc, n_inv);

        let mut a = Affine::<C>::zero();
        C::to_affine(&acc, &mut a);
        out.push(a);
    }

    out
}

pub fn construct_lagrange_basis<C, P>(n: usize, root: C::ScalarField) -> Vec<P>
where
    C: Curve,
    C::ScalarField: FieldImpl,
    <C as Curve>::ScalarField: Arithmetic,
    P: UnivariatePolynomial<Field = <C as Curve>::ScalarField>,
{
    let mut roots = Vec::with_capacity(n);
    let mut pow = C::ScalarField::one();
    for _ in 0..n {
        roots.push(pow);
        pow = pow * root;
    }
    
    let mut bases = Vec::with_capacity(n);
    
    for i in 0..n {
        let coeff = vec![C::ScalarField::one(), C::ScalarField::zero()];
        let mut li = P::from_coeffs(HostSlice::from_slice(&coeff), 2);

        let x_i = roots[i];

        for (j, &x_j) in roots.iter().enumerate() {
            if j == i { continue; }

            let bottom_inv = (x_i - x_j).inv();
            let c0 = C::ScalarField::zero() - (x_j * bottom_inv);
            let c1 = bottom_inv;
            let lii = P::from_coeffs(HostSlice::from_slice(&[c0, c1]), 2);
            li = li.mul(&lii);
        }

        bases.push(li);
    }
    bases
}

#[cfg(test)]
mod lagrange_test {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    use icicle_core::traits::FieldImpl;
    use icicle_core::ntt::{get_root_of_unity, initialize_domain, NTTInitDomainConfig, release_domain};
    use std::ops::Mul;

    use crate::{
        kzg::lagrange::{compute_lagrange_basis_commitments, construct_lagrange_basis},
        utils::srs::unsafe_setup_from_tau,
    };

    #[test]
    fn test_lagrange() {
        let n: usize = 16;
        let domain = get_root_of_unity::<Bn254ScalarField>((n * n * n * n).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        
        let root = get_root_of_unity::<Bn254ScalarField>((n).try_into().unwrap());
        
        let lb = construct_lagrange_basis::<Bn254CurveCfg, Bn254Poly>(n, root);

        let tau = Bn254ScalarField::from_u32(19u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(n - 1, tau);
        
        let mut lb_commitments_from_tau: Vec<Affine::<Bn254CurveCfg>> = Vec::with_capacity(lb.len());
        for i in 0..lb.len() {
            let li_tau = lb[i].eval(&tau);
            let projective_result = Bn254CurveCfg::get_generator().mul(li_tau);
            let mut affine_result = Affine::<Bn254CurveCfg>::zero();
            Bn254CurveCfg::to_affine(&projective_result, &mut affine_result);
            lb_commitments_from_tau.push(affine_result);
        }

        let srs_projective: Vec<Projective::<Bn254CurveCfg>> = srs.iter().map(|c| c.to_projective()).collect();
        let lb_commitments = compute_lagrange_basis_commitments(&srs_projective);
        
        release_domain::<Bn254ScalarField>().unwrap();

        assert_eq!(lb_commitments_from_tau, lb_commitments);
    }
}