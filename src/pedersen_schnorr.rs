use icicle_core::curve::{Curve, Affine};
use std::marker::PhantomData;
use icicle_core::traits::Arithmetic;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::GenerateRandom;

use self::structs::{Error, Instance, Proof, Witness};
use self::tr::Transcript;

pub mod structs;
mod tr;

/*
    Given points P and H prove knowledge of opening of pedersen commitments
    X = xQ + rH

    Round1:
    p.1.1. samples (b_1, b_2)
    p.1.2. sends R = b_1P + b_2H

    v.1.1 Sends random c

    Round2:
    p.2.1. sends z_1 = cx + b_1
    p.2.1. sends z_2 = cr + b_2

    v.2.1 cX + R = z_1P + z_2H
*/

pub struct Argument<C: Curve> {
    _c: PhantomData<C>,
}

impl<C: Curve> Argument<C> {
    pub fn prove (
        instance: &Instance<C>,
        witness: &Witness<C>,
    ) -> Proof<C> 
    where
        <C as Curve>::ScalarField: Arithmetic,
        <<C as Curve>::ScalarField as FieldImpl>::Config: GenerateRandom<<C as Curve>::ScalarField>
    {
        let mut tr = Transcript::<C>::new(b"pedersen-schnorr");
        tr.send_instance(instance);
        
        let random_scalars = <<C::ScalarField as FieldImpl>::Config as GenerateRandom<C::ScalarField>>::generate_random(2);
        let b_1 = random_scalars[0];
        let b_2 = random_scalars[1];

        let p_term = C::mul_scalar(instance.p_base.to_projective(), b_1);
        let h_term = C::mul_scalar(instance.h_base.to_projective(), b_2);
        let blinder = p_term + h_term;
        
        let mut blinder_affine = Affine::<C>::zero();
        C::to_affine(&blinder, &mut blinder_affine);

        tr.send_blinder(&blinder_affine);
        let c = tr.get_c();

        let z_1 = c * witness.x + b_1;
        let z_2 = c * witness.r + b_2;

        Proof { blinder: blinder_affine, z_1, z_2 }
    }

    pub fn verify(instance: &Instance<C>, proof: &Proof<C>) -> Result<(), Error> {
        let mut tr = Transcript::<C>::new(b"pedersen-schnorr");
        tr.send_instance(instance);

        tr.send_blinder(&proof.blinder);
        let c = tr.get_c();
        
        let x_projective = instance.x.to_projective();
        let blinder_projective = proof.blinder.to_projective();

        let x_c = C::mul_scalar(x_projective, c);
        let lhs = x_c + blinder_projective;

        let p_base_projective = instance.p_base.to_projective();
        let h_base_projective = instance.h_base.to_projective();

        let p_z_1 = C::mul_scalar(p_base_projective, proof.z_1);
        let h_z_2 = C::mul_scalar(h_base_projective, proof.z_2);
        let rhs = p_z_1 + h_z_2;

        if lhs != rhs {
            return Err(Error::RelationCheckFailed);
        }

        Ok(())
    }
}

#[cfg(test)]
mod pedersen_schnorr_test {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine};
    use icicle_core::traits::FieldImpl;

    use super::{
        structs::{Instance, Witness},
        Argument,
    };

    #[test]
    fn simple_test() {
        let g = Bn254CurveCfg::get_generator();

        // setup bases
        let p = Bn254ScalarField::from_u32(200u32);
        let h = Bn254ScalarField::from_u32(300u32);
        
        let p_base = Bn254CurveCfg::mul_scalar(g, p);
        let mut p_base_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&p_base, &mut p_base_affine);

        let h_base = Bn254CurveCfg::mul_scalar(g, h);
        let mut h_base_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&h_base, &mut h_base_affine);

        // witness
        let x_witness = Bn254ScalarField::from_u32(3u32);
        let r = Bn254ScalarField::from_u32(7u32);

        // instance
        let x = Bn254CurveCfg::mul_scalar(p_base, x_witness) + Bn254CurveCfg::mul_scalar(h_base, r);
        let mut x_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&x, &mut x_affine);

        let instance = Instance::<Bn254CurveCfg> {
            p_base: p_base_affine,
            h_base: h_base_affine,
            x: x_affine,
        };

        let witness = Witness::<Bn254CurveCfg> { x: x_witness, r };

        let proof = Argument::prove(&instance, &witness);
        let result = Argument::verify(&instance, &proof);
        assert!(result.is_ok());
    }
}