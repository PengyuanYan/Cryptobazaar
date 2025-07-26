use icicle_core::curve::{Curve,Affine};
use icicle_core::polynomials::UnivariatePolynomial;
use std::marker::PhantomData;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError, Valid};
use ark_std::io::{Read, Write};
use icicle_core::traits::FieldImpl;

use ark_serialize::Validate;
use ark_serialize::Compress;

#[derive(Debug)]
pub enum Error {
    OpeningFailed,
    RelationCheck,
}

pub struct Instance<C: Curve> {
    pub n: usize,
    pub mu: C::ScalarField,
    pub acc_cm: Affine::<C>,
}

// this is just for make the code can be compiled
impl<C: Curve> Valid for Instance<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.acc_cm.to_projective()) {
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
        writer.write_all(&self.mu.to_bytes_le())?;
        writer.write_all(&self.acc_cm.x.to_bytes_le())?;
        writer.write_all(&self.acc_cm.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len   = C::BaseField::zero().to_bytes_le().len();

        8 + scalar_len + 2 * base_len
    }
}

impl<C: Curve> CanonicalDeserialize for Instance<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len   = C::BaseField::zero().to_bytes_le().len();

        let mut n_buf = [0u8; 8];
        reader.read_exact(&mut n_buf)?;
        let n_u64 = u64::from_le_bytes(n_buf);
        let n_usize = usize::try_from(n_u64)
            .map_err(|_| SerializationError::InvalidData)?;

        let mut mu_bytes = vec![0u8; scalar_len];
        reader.read_exact(&mut mu_bytes)?;
        let mu = C::ScalarField::from_bytes_le(&mu_bytes);

        let mut x_bytes = vec![0u8; base_len];
        let mut y_bytes = vec![0u8; base_len];
        reader.read_exact(&mut x_bytes)?;
        reader.read_exact(&mut y_bytes)?;
        let x = C::BaseField::from_bytes_le(&x_bytes);
        let y = C::BaseField::from_bytes_le(&y_bytes);
        let acc_cm = Affine::<C> { x, y };

        if matches!(validate, Validate::Yes) && !C::is_on_curve(acc_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { n: n_usize, mu, acc_cm })
    }
}

pub struct Witness<C, P> 
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>,
{
    pub acc: P,
    pub _e: PhantomData<C>,
}

pub struct Proof<C: Curve> {
    pub q: C::ScalarField,
    pub acc_opening: C::ScalarField,
    pub acc_shifted_opening: C::ScalarField,

    pub q_0: Affine::<C>,
    pub q_1: Affine::<C>,
    pub q_2: Affine::<C>,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.q_0.to_projective()) {
            return Err(SerializationError::InvalidData);
        }
        if !C::is_on_curve(self.q_1.to_projective()) {
            return Err(SerializationError::InvalidData);
        }
        if !C::is_on_curve(self.q_2.to_projective()) {
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
        writer.write_all(&self.q.to_bytes_le())?;
        writer.write_all(&self.acc_opening.to_bytes_le())?;
        writer.write_all(&self.acc_shifted_opening.to_bytes_le())?;
        writer.write_all(&self.q_0.x.to_bytes_le())?;
        writer.write_all(&self.q_0.y.to_bytes_le())?;
        writer.write_all(&self.q_1.x.to_bytes_le())?;
        writer.write_all(&self.q_1.y.to_bytes_le())?;
        writer.write_all(&self.q_2.x.to_bytes_le())?;
        writer.write_all(&self.q_2.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len   = C::BaseField::zero().to_bytes_le().len();

        scalar_len * 3 + base_len * 6
    }
}

impl<C: Curve> CanonicalDeserialize for Proof<C> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress, // it cant change anything
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len   = C::BaseField::zero().to_bytes_le().len();
        
        let mut q_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_buf)?;
        let q = C::ScalarField::from_bytes_le(&q_buf);

        let mut opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut opening_buf )?;
        let opening = C::ScalarField::from_bytes_le(&opening_buf);

        let mut shifted_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut shifted_opening_buf)?;
        let shifted_opening = C::ScalarField::from_bytes_le(&shifted_opening_buf);

        let mut q_0_x_bytes = vec![0u8; base_len];
        let mut q_0_y_bytes = vec![0u8; base_len];
        reader.read_exact(&mut q_0_x_bytes)?;
        reader.read_exact(&mut q_0_y_bytes)?;
        let q_0_x = C::BaseField::from_bytes_le(&q_0_x_bytes);
        let q_0_y = C::BaseField::from_bytes_le(&q_0_y_bytes);
        let q_0 = Affine::<C> { x: q_0_x, y: q_0_y };
        
        let mut q_1_x_bytes = vec![0u8; base_len];
        let mut q_1_y_bytes = vec![0u8; base_len];
        reader.read_exact(&mut q_1_x_bytes)?;
        reader.read_exact(&mut q_1_y_bytes)?;
        let q_1_x = C::BaseField::from_bytes_le(&q_1_x_bytes);
        let q_1_y = C::BaseField::from_bytes_le(&q_1_y_bytes);
        let q_1 = Affine::<C> { x: q_1_x, y: q_1_y };

        let mut q_2_x_bytes = vec![0u8; base_len];
        let mut q_2_y_bytes = vec![0u8; base_len];
        reader.read_exact(&mut q_2_x_bytes)?;
        reader.read_exact(&mut q_2_y_bytes)?;
        let q_2_x = C::BaseField::from_bytes_le(&q_2_x_bytes);
        let q_2_y = C::BaseField::from_bytes_le(&q_2_y_bytes);
        let q_2 = Affine::<C> { x: q_2_x, y: q_2_y };

        if matches!(validate, Validate::Yes) 
            && !C::is_on_curve(q_0.to_projective())
            && !C::is_on_curve(q_1.to_projective())
            && !C::is_on_curve(q_2.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { q: q, acc_opening: opening, acc_shifted_opening: shifted_opening, q_0: q_0, q_1: q_1, q_2: q_2 })
    }
}