use icicle_core::traits::FieldImpl;
use icicle_core::ntt::NTTDomain;
use icicle_core::curve::{Curve,Affine};
use icicle_core::polynomials::UnivariatePolynomial;
use std::marker::PhantomData;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;
use ark_serialize::Valid;
use ark_serialize::SerializationError;


#[derive(Debug)]
pub enum Error {
    WellFormation,
    Sumcheck,
    Openings,
}

pub struct VerifierIndex<C: Curve> {
    pub s_cm: Affine::<C>,
}

impl<C: Curve> Valid for VerifierIndex<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.s_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }
        Ok(())
    }
}

impl<C: Curve> CanonicalSerialize for VerifierIndex<C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&self.s_cm.x.to_bytes_le())?;
        writer.write_all(&self.s_cm.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        
        2 * base_len
    }
}

impl<C: Curve> CanonicalDeserialize for VerifierIndex<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let base_len = C::BaseField::zero().to_bytes_le().len();

        let mut x_bytes = vec![0u8; base_len];
        let mut y_bytes = vec![0u8; base_len];
        reader.read_exact(&mut x_bytes)?;
        reader.read_exact(&mut y_bytes)?;
        let x = C::BaseField::from_bytes_le(&x_bytes);
        let y = C::BaseField::from_bytes_le(&y_bytes);
        let s_cm = Affine::<C> { x, y };

        if matches!(validate, Validate::Yes) && !C::is_on_curve(s_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { s_cm })
    }
}

pub struct ProverIndex<C, P>
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>
{
    pub s: P,
    pub s_coset_evals: Vec<C::ScalarField>,
}

pub struct Instance<C: Curve> {
    pub f_cm: Affine::<C>,
}

impl<C: Curve> Valid for Instance<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.f_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }
        Ok(())
    }
}

impl<C: Curve> CanonicalSerialize for Instance<C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&self.f_cm.x.to_bytes_le())?;
        writer.write_all(&self.f_cm.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        
        2 * base_len
    }
}

impl<C: Curve> CanonicalDeserialize for Instance<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let base_len = C::BaseField::zero().to_bytes_le().len();

        let mut x_bytes = vec![0u8; base_len];
        let mut y_bytes = vec![0u8; base_len];
        reader.read_exact(&mut x_bytes)?;
        reader.read_exact(&mut y_bytes)?;
        let x = C::BaseField::from_bytes_le(&x_bytes);
        let y = C::BaseField::from_bytes_le(&y_bytes);
        let f_cm = Affine::<C> { x, y };

        if matches!(validate, Validate::Yes) && !C::is_on_curve(f_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { f_cm })
    }
}

pub struct Witness<C, P>
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>,
{
    pub f: P,
    pub e: PhantomData<C>,
}

pub struct Proof<C: Curve> {
    pub gamma: C::ScalarField,

    pub b_cm: Affine::<C>,
    pub q_cm: Affine::<C>,

    pub f_opening: C::ScalarField,
    pub s_opening: C::ScalarField,
    pub b_opening: C::ScalarField,
    pub q_opening: C::ScalarField,

    pub q_0: Affine::<C>,
    pub q_1: Affine::<C>,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.b_cm.to_projective()) &&
           !C::is_on_curve(self.q_cm.to_projective()) &&

           !C::is_on_curve(self.q_0.to_projective()) &&
           !C::is_on_curve(self.q_1.to_projective())
        
        {
            return Err(SerializationError::InvalidData);
        }
        Ok(())
    }
}

impl<C: Curve> CanonicalSerialize for Proof<C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&self.gamma.to_bytes_le())?;

        writer.write_all(&self.b_cm.x.to_bytes_le())?;
        writer.write_all(&self.b_cm.y.to_bytes_le())?;

        writer.write_all(&self.q_cm.x.to_bytes_le())?;
        writer.write_all(&self.q_cm.y.to_bytes_le())?;

        writer.write_all(&self.f_opening.to_bytes_le())?;
        writer.write_all(&self.s_opening.to_bytes_le())?;
        writer.write_all(&self.b_opening.to_bytes_le())?;
        writer.write_all(&self.q_opening.to_bytes_le())?;
        
        writer.write_all(&self.q_0.x.to_bytes_le())?;
        writer.write_all(&self.q_0.y.to_bytes_le())?;
        writer.write_all(&self.q_1.x.to_bytes_le())?;
        writer.write_all(&self.q_1.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        5 * scalar_len + 2 * base_len * 4
    }
}

impl<C: Curve> CanonicalDeserialize for Proof<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();
        
        let mut gamma_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut gamma_buf)?;
        let gamma = C::ScalarField::from_bytes_le(&gamma_buf);

        let mut read_affine = |reader: &mut dyn Read| -> Result<Affine::<C>, std::io::Error> {
            let mut x_bytes = vec![0u8; base_len];
            let mut y_bytes = vec![0u8; base_len];
            reader.read_exact(&mut x_bytes)?;
            reader.read_exact(&mut y_bytes)?;
            let x = C::BaseField::from_bytes_le(&x_bytes);
            let y = C::BaseField::from_bytes_le(&y_bytes);
            Ok(Affine::<C> { x, y })
        };

        let b_cm = read_affine(&mut reader)?;
        let q_cm = read_affine(&mut reader)?;

        let mut f_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut f_opening_buf)?;
        let f_opening = C::ScalarField::from_bytes_le(&f_opening_buf);

        let mut s_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut s_opening_buf)?;
        let s_opening = C::ScalarField::from_bytes_le(&s_opening_buf);

        let mut b_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut b_opening_buf)?;
        let b_opening = C::ScalarField::from_bytes_le(&b_opening_buf);

        let mut q_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_opening_buf)?;
        let q_opening = C::ScalarField::from_bytes_le(&q_opening_buf);

        let q_0 = read_affine(&mut reader)?;
        let q_1 = read_affine(&mut reader)?;

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(b_cm.to_projective()) &&
            !C::is_on_curve(q_cm.to_projective()) &&
            !C::is_on_curve(q_0.to_projective()) &&
            !C::is_on_curve(q_1.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { gamma,
                  b_cm, q_cm,
                  f_opening, s_opening, b_opening, q_opening,
                  q_0, q_1,
        })
    }
}

#[cfg(test)]
mod serialize_test {
    use std::ops::Mul;
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::traits::FieldImpl;
    use rand_chacha::ChaCha20Rng;
    use std::marker::PhantomData;
    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    
    use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, Compress};
    use ark_serialize::Validate;

    use crate::bid_encoder::BidEncoder;
    use crate::{
        kzg::{PK, VK},
        utils::srs::unsafe_setup_from_tau,
    };

    use super::{Proof};

    #[test]
    fn test_serialize() {
        let gamma = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let b_cm = Affine::<Bn254CurveCfg>::zero();
        let q_cm = Affine::<Bn254CurveCfg>::zero();
        
        let f_opening = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let s_opening = Bn254ScalarField::one();
        let b_opening = Bn254ScalarField::one() + Bn254ScalarField::one();
        let q_opening = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();

        let q_0 = Affine::<Bn254CurveCfg>::zero();
        let q_1 = Affine::<Bn254CurveCfg>::zero();

        let proof = Proof { gamma,
                            b_cm, q_cm,
                            f_opening, s_opening, b_opening, q_opening,
                            q_0, q_1,
                           };
        
        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No);

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.gamma, gamma);
        assert_eq!(result.b_cm, b_cm);
        assert_eq!(result.b_opening, b_opening);
    }
}