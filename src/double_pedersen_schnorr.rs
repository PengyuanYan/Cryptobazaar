use icicle_core::curve::{Curve, Projective, Affine};
use std::marker::PhantomData;
use icicle_core::traits::Arithmetic;
use icicle_core::traits::FieldImpl;
use icicle_core::traits::GenerateRandom;

use self::structs::{Error, Instance, Proof, Witness};
use self::tr::Transcript;

pub mod structs;
mod tr;

/*
    Given points Q, P and H prove knowledge of opening of pedersen commitments
    1. X_1 = aQ + r_1H
    2. X_2 = aP + r_2H

    Round1:
    p.1.1. samples (b_1, b_2, b_3)
    p.1.2. sends R_1 = b_1Q + b_2H
    p.1.3. sends R_2 = b_1P + b_3H

    v.1.1 Sends random c

    Round2:
    p.2.1. sends z_1 = ca + b_1
    p.2.1. sends z_2 = cr_1 + b_2
    p.2.1. sends z_3 = cr_2 + b_3

    v.2.1 cX_1 + R_1 = z_1Q + z_2H
    v.2.2 cX_2 + R_2 = z_1P + z_3H
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

        let (b_1, b_2, b_3) = (
            <<C::ScalarField as FieldImpl>::Config as GenerateRandom<C::ScalarField>>::generate_random(1)[0],
            <<C::ScalarField as FieldImpl>::Config as GenerateRandom<C::ScalarField>>::generate_random(1)[0],
            <<C::ScalarField as FieldImpl>::Config as GenerateRandom<C::ScalarField>>::generate_random(1)[0],
        );

        let rand_1 = C::mul_scalar(instance.q_base.to_projective(), b_1) + C::mul_scalar(instance.h_base.to_projective(), b_2);
        let rand_2 = C::mul_scalar(instance.p_base.to_projective(), b_1) + C::mul_scalar(instance.h_base.to_projective(), b_3);
        
        let mut rand_1_affine = Affine::<C>::zero();
        C::to_affine(&rand_1, &mut rand_1_affine);
        let mut rand_2_affine = Affine::<C>::zero();
        C::to_affine(&rand_2, &mut rand_2_affine);

        tr.send_blinders(&rand_1_affine, &rand_2_affine);
        let c = tr.get_c();

        let z_1 = c * witness.a + b_1;
        let z_2 = c * witness.r_1 + b_2;
        let z_3 = c * witness.r_2 + b_3;

        Proof {
            rand_1: rand_1_affine,
            rand_2: rand_2_affine,
            z_1,
            z_2,
            z_3,
        }
    }

    pub fn verify(instance: &Instance<C>, proof: &Proof<C>) -> Result<(), Error> {
        let mut tr = Transcript::<C>::new(b"pedersen-schnorr");
        tr.send_instance(instance);

        tr.send_blinders(&proof.rand_1, &proof.rand_2);
        let c = tr.get_c();
        
        let x_1_projective = instance.x_1.to_projective();
        let rand_1_projective = proof.rand_1.to_projective();

        let x_1_c = C::mul_scalar(x_1_projective, c);
        let lhs_1 = x_1_c + rand_1_projective;

        let q_base_projective = instance.q_base.to_projective();
        let h_base_projective = instance.h_base.to_projective();

        let q_z_1 = C::mul_scalar(q_base_projective, proof.z_1);
        let h_z_2 = C::mul_scalar(h_base_projective, proof.z_2);
        let rhs_1 = q_z_1 + h_z_2;

        let eq1 = { lhs_1 == rhs_1 };

        if !eq1 {
            return Err(Error::Relation1);
        }
        
        let x_2_projective = instance.x_2.to_projective();
        let rand_2_projective = proof.rand_2.to_projective();

        let x_2_c = C::mul_scalar(x_2_projective, c);
        let lhs_2 = x_2_c + rand_2_projective;

        let p_base_projective = instance.p_base.to_projective();
        let h_base_projective = instance.h_base.to_projective();

        let p_z_1 = C::mul_scalar(p_base_projective, proof.z_1);
        let h_z_3 = C::mul_scalar(h_base_projective, proof.z_3);
        let rhs_2 = p_z_1 + h_z_3;

        let eq2 = { lhs_2 == rhs_2 };

        if !eq2 {
            return Err(Error::Relation1);
        }

        Ok(())
    }
}

#[cfg(test)]
mod pedersen_schnorr_test {
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::traits::FieldImpl;
    use icicle_runtime::memory::HostSlice;
    use icicle_core::traits::GenerateRandom;

    use super::{
        structs::{Instance, Witness},
        Argument,
    };

    #[test]
    fn simple_test() {
        let g = Bn254CurveCfg::get_generator();

        // setup bases
        let q = Bn254ScalarField::from_u32(100u32);
        let p = Bn254ScalarField::from_u32(200u32);
        let h = Bn254ScalarField::from_u32(300u32);

        let q_base = Bn254CurveCfg::mul_scalar(g, q);
        let p_base = Bn254CurveCfg::mul_scalar(g, p);
        let h_base = Bn254CurveCfg::mul_scalar(g, h);

        // witness
        let a = Bn254ScalarField::from_u32(3u32);
        let r_1 = Bn254ScalarField::from_u32(7u32);
        let r_2 = Bn254ScalarField::from_u32(13u32);

        // instance
        let x_1 = Bn254CurveCfg::mul_scalar(q_base, a) + Bn254CurveCfg::mul_scalar(h_base, r_1);
        let x_2 = Bn254CurveCfg::mul_scalar(p_base, a) + Bn254CurveCfg::mul_scalar(h_base, r_2);

        let instance = Instance::<Bn254CurveCfg> {
            q_base: q_base.into(),
            p_base: p_base.into(),
            h_base: h_base.into(),
            x_1: x_1.into(),
            x_2: x_2.into(),
        };

        let witness = Witness::<Bn254CurveCfg> { a, r_1, r_2 };

        let proof = Argument::prove(&instance, &witness);
        let result = Argument::verify(&instance, &proof);
        assert!(result.is_ok());
    }
}