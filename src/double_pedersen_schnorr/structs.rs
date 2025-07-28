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
    Relation1,
    Relation2,
}

pub struct Instance<C: Curve> {
    pub q_base: Affine::<C>,
    pub p_base: Affine::<C>,
    pub h_base: Affine::<C>,

    pub x_1: Affine::<C>,
    pub x_2: Affine::<C>,
}

impl<C: Curve> Valid for Instance<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.q_base.to_projective()) &&
           !C::is_on_curve(self.p_base.to_projective()) &&
           !C::is_on_curve(self.h_base.to_projective()) &&
           !C::is_on_curve(self.x_1.to_projective()) &&
           !C::is_on_curve(self.x_2.to_projective())
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
        writer.write_all(&self.q_base.x.to_bytes_le())?;
        writer.write_all(&self.q_base.y.to_bytes_le())?;

        writer.write_all(&self.p_base.x.to_bytes_le())?;
        writer.write_all(&self.p_base.y.to_bytes_le())?;

        writer.write_all(&self.h_base.x.to_bytes_le())?;
        writer.write_all(&self.h_base.y.to_bytes_le())?;

        writer.write_all(&self.x_1.x.to_bytes_le())?;
        writer.write_all(&self.x_1.y.to_bytes_le())?;

        writer.write_all(&self.x_2.x.to_bytes_le())?;
        writer.write_all(&self.x_2.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();

        2 * base_len * 5
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
        
        let q_base = read_affine(&mut reader)?;
        let p_base = read_affine(&mut reader)?;
        let h_base = read_affine(&mut reader)?;
        let x_1 = read_affine(&mut reader)?;
        let x_2 = read_affine(&mut reader)?; 

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(q_base.to_projective()) &&
            !C::is_on_curve(p_base.to_projective()) &&
            !C::is_on_curve(h_base.to_projective()) &&
            !C::is_on_curve(x_1.to_projective()) &&
            !C::is_on_curve(x_2.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { q_base, p_base, h_base,
                  x_1, x_2,
        })
    }
}

pub struct Witness<C: Curve> {
    pub a: C::ScalarField,
    pub r_1: C::ScalarField,
    pub r_2: C::ScalarField,
}

pub struct Proof<C: Curve> {
    pub rand_1: Affine::<C>,
    pub rand_2: Affine::<C>,

    pub z_1: C::ScalarField,
    pub z_2: C::ScalarField,
    pub z_3: C::ScalarField,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.rand_1.to_projective()) &&
           !C::is_on_curve(self.rand_2.to_projective())
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
        writer.write_all(&self.rand_1.x.to_bytes_le())?;
        writer.write_all(&self.rand_1.y.to_bytes_le())?;
        writer.write_all(&self.rand_2.x.to_bytes_le())?;
        writer.write_all(&self.rand_2.y.to_bytes_le())?;

        writer.write_all(&self.z_1.to_bytes_le())?;
        writer.write_all(&self.z_2.to_bytes_le())?;
        writer.write_all(&self.z_3.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        3 * scalar_len + 2 * base_len * 2
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
        
        let rand_1 = read_affine(&mut reader)?;
        let rand_2 = read_affine(&mut reader)?;

        let mut z_1_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_1_buf)?;
        let z_1 = C::ScalarField::from_bytes_le(&z_1_buf);
        let mut z_2_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_2_buf)?;
        let z_2 = C::ScalarField::from_bytes_le(&z_2_buf);
        let mut z_3_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_3_buf)?;
        let z_3 = C::ScalarField::from_bytes_le(&z_3_buf);

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(rand_1.to_projective()) &&
            !C::is_on_curve(rand_2.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { rand_1, rand_2,
                  z_1, z_2, z_3,
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
        let q_base = Affine::<Bn254CurveCfg>::zero();
        let p_base = Affine::<Bn254CurveCfg>::zero();
        let h_base = Affine::<Bn254CurveCfg>::zero();
        let x_1 = Affine::<Bn254CurveCfg>::zero();
        let x_2 = Affine::<Bn254CurveCfg>::zero();

        let instance = Instance { q_base, p_base, h_base,
                                  x_1, x_2,
                                };
        
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Instance::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(instance.q_base, result.q_base);
        assert_eq!(instance.p_base, result.p_base);
        assert_eq!(instance.h_base, result.h_base);
        assert_eq!(instance.x_1, result.x_1);
        assert_eq!(instance.x_2, result.x_2);

        let rand_1 = Affine::<Bn254CurveCfg>::zero();
        let rand_2 = Affine::<Bn254CurveCfg>::zero();

        let z_1 = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let z_2 = Bn254ScalarField::one();
        let z_3 = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();

        let proof = Proof { rand_1, rand_2,
                            z_1, z_2, z_3,
                           };
        
        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.rand_1, proof.rand_1);
        assert_eq!(result.z_1, proof.z_1);
        assert_eq!(result.z_2, proof.z_2);
        assert_eq!(result.z_3, proof.z_3);
    }
}