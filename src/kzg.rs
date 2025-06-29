use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use std::marker::PhantomData;

use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::HostSlice;

use icicle_core::polynomials::UnivariatePolynomial;
use crate::utils::get_coeffs_of_poly;

pub struct Kzg<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    _e: PhantomData<(C1, C2, F)>,
}

pub struct PK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub srs: Vec<Affine<C1>>,
    _e: PhantomData<(C2, F)>,
}

pub struct VK<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
{
    pub g2: C2,
    pub neg_x_g2: C2,
    _e: PhantomData<(C1, F)>,
}

impl<C1, C2, F> Kzg<C1, C2, F>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    C1: icicle_core::msm::MSM<C1>,
{
    pub fn commit<P>(
        pk: &PK<C1, C2, F>,
        poly: &P
    ) -> Affine<C1> 
    where
        P: UnivariatePolynomial<Field = C1::ScalarField>,
    {
        if pk.srs.len() - 1 < poly.degree().try_into().unwrap() {
            panic!(
                "SRS size too small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                pk.srs.len()
            );
        }

        let mut cfg = MSMConfig::default();
        let mut projective_output = vec![Projective::<C1>::zero(); 1];
        let coeffs = get_coeffs_of_poly(poly);

        msm::msm(
            HostSlice::from_slice(&coeffs),
            HostSlice::from_slice(&pk.srs),
            &cfg,
            HostSlice::from_mut_slice(&mut projective_output),
        )
        .unwrap();
        
        let mut affine_output = Affine::<C1>::zero();
        C1::to_affine(&projective_output[0], &mut affine_output);
        affine_output
    }
}

#[cfg(test)]
mod test_kzg {
    use icicle_bn254::polynomials::DensePolynomial as Bn254DensePolynomial;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;

    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_runtime::memory::HostSlice;
    
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg, G1Projective, G1Affine};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_core::curve::Curve;

    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;
    
    use icicle_core::{msm, msm::MSMConfig};
    
    use crate::utils::srs::unsafe_setup_from_tau;

    use super::{Kzg, PK, VK};
    use std::marker::PhantomData;
    
    #[test]
    fn test_kzg_real() {
        let size = 16;
        let coeffs = ScalarCfg::generate_random(size);
        let poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&coeffs), size);
        let tau = Bn254ScalarField::from_u32(100u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(size - 1, tau);
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), _e: PhantomData, };

        let mut cfg = MSMConfig::default();
        let mut output = vec![G1Projective::zero(); size];
        
        let commit = Kzg::commit(&pk, &poly);
    }

    #[test]
    fn test_kzg_related_functions() {
        let size = 16;
        let coeffs = ScalarCfg::generate_random(size);
        let poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&coeffs), size);
        let tau = Bn254ScalarField::from_u32(100u32);
        println!("{tau}\n");

        let mut scalars = Vec::with_capacity(size);
        let mut acc = Bn254ScalarField::from_u32(1u32);
        for _ in 0..size {
            scalars.push(acc);
            acc = acc * tau;
            println!("{acc}");
        }
        assert_eq!(scalars.len(),size);

        println!("{scalars:?}");
        
        let projective_g1: G1Projective = Bn254CurveCfg::get_generator();

        println!("\n{projective_g1:?}");
        
        let mut srs = Vec::with_capacity(size);
        for i in &scalars {
            let projective_base = Bn254CurveCfg::mul_scalar(projective_g1, *i);
            let mut affine_base: G1Affine = G1Affine::zero();
            Bn254CurveCfg::to_affine(&projective_base, &mut affine_base);
            srs.push(affine_base);
        }
        
        println!("\n\n\n\n\n\n\n\n\n{srs:?}\n\n\n\n\n\n\n\n\n");
        
        let mut cfg = MSMConfig::default();
        let mut output = vec![G1Projective::zero(); size];
        
        msm::msm(
            HostSlice::from_slice(&coeffs),
            HostSlice::from_slice(&srs),
            &cfg,
            HostSlice::from_mut_slice(&mut output),
        )
        .unwrap();

        println!("\n{output:?}");

        let w = Bn254ScalarField::from_u32(10u32);

        let eval = poly.eval(&w);
    } 
}