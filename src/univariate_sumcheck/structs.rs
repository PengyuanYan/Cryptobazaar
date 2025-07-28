use icicle_core::traits::FieldImpl;
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
    Error,
}

/// An instance for the univariate sumcheck argument,
/// for some witness a, b such that:
/// - a_cm = Comm(a),
/// - b_cm = Comm(b),
/// - sum = ∑ a_i • b_i
pub struct Instance<C: Curve> {
    pub(crate) n: usize,
    pub(crate) a_cm: Affine::<C>,
    pub(crate) b_cm: Affine::<C>,
    pub(crate) sum: C::ScalarField,
}

impl<C: Curve> Valid for Instance<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.a_cm.to_projective()) &&
           !C::is_on_curve(self.b_cm.to_projective())        
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
        writer.write_all(&(self.n as u64).to_le_bytes())?;

        writer.write_all(&self.a_cm.x.to_bytes_le())?;
        writer.write_all(&self.a_cm.y.to_bytes_le())?;

        writer.write_all(&self.b_cm.x.to_bytes_le())?;
        writer.write_all(&self.b_cm.y.to_bytes_le())?;

        writer.write_all(&self.sum.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        8 + scalar_len + 2 * base_len * 2
    }
}

impl<C: Curve> CanonicalDeserialize for Instance<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();
        
        let mut n_buf = [0u8; 8];
        reader.read_exact(&mut n_buf)?;
        let n_u64 = u64::from_le_bytes(n_buf);
        let n_usize = usize::try_from(n_u64)
            .map_err(|_| SerializationError::InvalidData)?;
        
        let read_affine = |reader: &mut dyn Read| -> Result<Affine::<C>, std::io::Error> {
            let mut x_bytes = vec![0u8; base_len];
            let mut y_bytes = vec![0u8; base_len];
            reader.read_exact(&mut x_bytes)?;
            reader.read_exact(&mut y_bytes)?;
            let x = C::BaseField::from_bytes_le(&x_bytes);
            let y = C::BaseField::from_bytes_le(&y_bytes);
            Ok(Affine::<C> { x, y })
        };
        
        let a_cm = read_affine(&mut reader)?;
        let b_cm = read_affine(&mut reader)?;

        let mut sum_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut sum_buf)?;
        let sum = C::ScalarField::from_bytes_le(&sum_buf);


        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(a_cm.to_projective()) &&
            !C::is_on_curve(b_cm.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { n: n_usize,
                  a_cm, b_cm,
                  sum,
        })
    }
}

pub struct Witness<C, P> 
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>,
{
    pub(crate) a_poly: P,
    pub(crate) b_poly: P,
    pub e: PhantomData<C>,
}

pub struct Proof<C: Curve> {
    pub(crate) r_cm: Affine::<C>,
    pub(crate) q_cm: Affine::<C>,
    // opening proofs
    pub(crate) a_opening: C::ScalarField,
    pub(crate) b_opening: C::ScalarField,
    pub(crate) r_opening: C::ScalarField,
    pub(crate) q_opening: C::ScalarField,

    pub(crate) batch_opening_proof: Affine::<C>,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.r_cm.to_projective()) &&
           !C::is_on_curve(self.q_cm.to_projective()) &&
           !C::is_on_curve(self.batch_opening_proof.to_projective())
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
        writer.write_all(&self.r_cm.x.to_bytes_le())?;
        writer.write_all(&self.r_cm.y.to_bytes_le())?;
        writer.write_all(&self.q_cm.x.to_bytes_le())?;
        writer.write_all(&self.q_cm.y.to_bytes_le())?;

        writer.write_all(&self.a_opening.to_bytes_le())?;
        writer.write_all(&self.b_opening.to_bytes_le())?;
        writer.write_all(&self.r_opening.to_bytes_le())?;
        writer.write_all(&self.q_opening.to_bytes_le())?;

        writer.write_all(&self.batch_opening_proof.x.to_bytes_le())?;
        writer.write_all(&self.batch_opening_proof.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        4 * scalar_len + 2 * base_len * 3
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
        
        let r_cm = read_affine(&mut reader)?;
        let q_cm = read_affine(&mut reader)?;

        let mut a_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut a_opening_buf)?;
        let a_opening = C::ScalarField::from_bytes_le(&a_opening_buf);
        let mut b_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut b_opening_buf)?;
        let b_opening = C::ScalarField::from_bytes_le(&b_opening_buf);
        let mut r_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut r_opening_buf)?;
        let r_opening = C::ScalarField::from_bytes_le(&r_opening_buf);
        let mut q_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_opening_buf)?;
        let q_opening = C::ScalarField::from_bytes_le(&q_opening_buf);

        let batch_opening_proof = read_affine(&mut reader)?;

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(r_cm.to_projective()) &&
            !C::is_on_curve(q_cm.to_projective()) &&
            !C::is_on_curve(batch_opening_proof.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { r_cm, q_cm,
                  a_opening, b_opening, r_opening, q_opening,
                  batch_opening_proof,
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
        let n = 7usize;
        let a_cm = Affine::<Bn254CurveCfg>::zero();
        let b_cm = Affine::<Bn254CurveCfg>::zero();
        let sum = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();

        let instance = Instance { n,
                                  a_cm, b_cm,
                                  sum,
                                };
        
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Instance::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(instance.n, result.n);
        assert_eq!(instance.a_cm, result.a_cm);
        assert_eq!(instance.b_cm, result.b_cm);
        assert_eq!(instance.sum, result.sum);

        let r_cm = Affine::<Bn254CurveCfg>::zero();
        let q_cm = Affine::<Bn254CurveCfg>::zero();

        let a_opening = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let b_opening = Bn254ScalarField::one();
        let r_opening = Bn254ScalarField::one() + Bn254ScalarField::one();
        let q_opening = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();

        let batch_opening_proof = Affine::<Bn254CurveCfg>::zero();

        let proof = Proof { r_cm, q_cm,
                            a_opening, b_opening, r_opening, q_opening,
                            batch_opening_proof,
                           };
        
        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.a_opening, proof.a_opening);
        assert_eq!(result.b_opening, proof.b_opening);
        assert_eq!(result.r_opening, proof.r_opening);
        assert_eq!(result.q_opening, proof.q_opening);
    }
}