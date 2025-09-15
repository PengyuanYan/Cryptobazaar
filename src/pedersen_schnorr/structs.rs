// This file contains the Structs for pedersen_schnorr and corresponding serlization functions.
use icicle_core::traits::FieldImpl;
use icicle_core::curve::{Curve,Affine};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;
use ark_serialize::Valid;
use ark_serialize::SerializationError;

#[derive(Debug)]
pub enum Error {
    RelationCheckFailed,
}

pub struct Instance<C: Curve> {
    pub p_base: Affine::<C>,
    pub h_base: Affine::<C>,

    pub x: Affine::<C>,
}

impl<C: Curve> Valid for Instance<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.p_base.to_projective()) &&
           !C::is_on_curve(self.h_base.to_projective()) &&
           !C::is_on_curve(self.x.to_projective()) 
        {
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
        writer.write_all(&self.p_base.x.to_bytes_le())?;
        writer.write_all(&self.p_base.y.to_bytes_le())?;

        writer.write_all(&self.h_base.x.to_bytes_le())?;
        writer.write_all(&self.h_base.y.to_bytes_le())?;

        writer.write_all(&self.x.x.to_bytes_le())?;
        writer.write_all(&self.x.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();

        2 * base_len * 3
    }
}

impl<C: Curve> CanonicalDeserialize for Instance<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        
        let read_affine = |reader: &mut dyn Read| -> Result<Affine::<C>, std::io::Error> {
            let mut x_bytes = vec![0u8; base_len];
            let mut y_bytes = vec![0u8; base_len];
            reader.read_exact(&mut x_bytes)?;
            reader.read_exact(&mut y_bytes)?;
            let x = C::BaseField::from_bytes_le(&x_bytes);
            let y = C::BaseField::from_bytes_le(&y_bytes);
            Ok(Affine::<C> { x, y })
        };
        
        let p_base = read_affine(&mut reader)?;
        let h_base = read_affine(&mut reader)?;
        let x = read_affine(&mut reader)?; 

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(p_base.to_projective()) &&
            !C::is_on_curve(h_base.to_projective()) &&
            !C::is_on_curve(x.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { p_base,
                  h_base,
                  x,
        })
    }
}

pub struct Witness<C: Curve> {
    pub x: C::ScalarField,
    pub r: C::ScalarField,
}

pub struct Proof<C: Curve> {
    pub blinder: Affine::<C>,

    pub z_1: C::ScalarField,
    pub z_2: C::ScalarField,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.blinder.to_projective())
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
        writer.write_all(&self.blinder.x.to_bytes_le())?;
        writer.write_all(&self.blinder.y.to_bytes_le())?;

        writer.write_all(&self.z_1.to_bytes_le())?;
        writer.write_all(&self.z_2.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        2 * scalar_len + 2 * base_len
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
        
        let read_affine = |reader: &mut dyn Read| -> Result<Affine::<C>, std::io::Error> {
            let mut x_bytes = vec![0u8; base_len];
            let mut y_bytes = vec![0u8; base_len];
            reader.read_exact(&mut x_bytes)?;
            reader.read_exact(&mut y_bytes)?;
            let x = C::BaseField::from_bytes_le(&x_bytes);
            let y = C::BaseField::from_bytes_le(&y_bytes);
            Ok(Affine::<C> { x, y })
        };
        
        let blinder = read_affine(&mut reader)?;

        let mut z_1_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_1_buf)?;
        let z_1 = C::ScalarField::from_bytes_le(&z_1_buf);
        let mut z_2_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_2_buf)?;
        let z_2 = C::ScalarField::from_bytes_le(&z_2_buf);

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(blinder.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { blinder,
                  z_1, z_2
        })
    }
}

#[cfg(test)]
mod serialize_test {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::Affine;
    use icicle_core::traits::FieldImpl;
    
    use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, Compress};
    use ark_serialize::Validate;

    use super::{Instance, Proof};

    #[test]
    fn test_serialize() {
        let p_base = Affine::<Bn254CurveCfg>::zero();
        let h_base = Affine::<Bn254CurveCfg>::zero();
        let x = Affine::<Bn254CurveCfg>::zero();

        let instance = Instance { p_base,
                                  h_base,
                                  x,
                                };
        
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Instance::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(instance.p_base, result.p_base);
        assert_eq!(instance.h_base, result.h_base);
        assert_eq!(instance.x, result.x);

        let blinder = Affine::<Bn254CurveCfg>::zero();

        let z_1 = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let z_2 = Bn254ScalarField::one();

        let proof = Proof { blinder,
                            z_1, z_2,
                           };
        
        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.blinder, proof.blinder);
        assert_eq!(result.z_1, proof.z_1);
        assert_eq!(result.z_2, proof.z_2);
    }
}