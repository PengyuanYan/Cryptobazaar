use icicle_core::traits::FieldImpl;
use icicle_core::curve::{Curve,Affine};
use icicle_core::polynomials::UnivariatePolynomial;
use std::marker::PhantomData;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError, Valid};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;

#[derive(Debug)]
pub enum Error {
    Opening,
    ShiftedOpening,
    RelationCheck,
}

pub struct Oracle<'a, F: FieldImpl>(pub(crate) &'a [F]);

pub struct VerifierIndex<C: Curve> {
    pub q_price_cm: Affine::<C>,
}

// this is just for make the code can be compiled
impl<C: Curve> Valid for VerifierIndex<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.q_price_cm.to_projective()) {
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
        writer.write_all(&self.q_price_cm.x.to_bytes_le())?;
        writer.write_all(&self.q_price_cm.y.to_bytes_le())?;

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
        let q_price_cm = Affine::<C> { x, y };

        if matches!(validate, Validate::Yes) && !C::is_on_curve(q_price_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { q_price_cm })
    }
}

pub struct ProverIndex<C, P>
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>,
{
    pub q_price: P,
    pub q_price_coset_evals: Vec<C::ScalarField>,
    pub l_p_coset_evals: Vec<C::ScalarField>,
}

pub struct Witness<C, P>
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>,
{
    pub bid: P,
    pub f: P,
    pub r: P,
    pub r_inv: P,
    pub diff: P,
    pub g: P,
    pub e: PhantomData<C>,
}

//#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<C: Curve> {
    pub bid_cm: Affine::<C>,
    pub r_cm: Affine::<C>,
    pub r_inv_cm: Affine::<C>,
    pub f_cm: Affine::<C>,
    pub diff_cm: Affine::<C>,
    pub g_cm: Affine::<C>,

    pub q_price_opening: C::ScalarField,
    pub bid_opening: C::ScalarField,
    pub bid_shift_opening: C::ScalarField,
    pub f_opening: C::ScalarField,
    pub r_opening: C::ScalarField,
    pub r_inv_opening: C::ScalarField,
    pub diff_opening: C::ScalarField,
    pub g_opening: C::ScalarField,

    pub q_chunk_0_cm: Affine::<C>,
    pub q_chunk_1_cm: Affine::<C>,

    pub q_chunk_0_opening: C::ScalarField,
    pub q_chunk_1_opening: C::ScalarField,

    pub w_0: Affine::<C>,
    pub w_1: Affine::<C>,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.bid_cm.to_projective()) &&
           !C::is_on_curve(self.r_cm.to_projective()) &&
           !C::is_on_curve(self.r_inv_cm.to_projective()) &&
           !C::is_on_curve(self.f_cm.to_projective()) &&
           !C::is_on_curve(self.diff_cm.to_projective()) &&
           !C::is_on_curve(self.g_cm.to_projective()) &&

           !C::is_on_curve(self.q_chunk_0_cm.to_projective()) &&
           !C::is_on_curve(self.q_chunk_1_cm.to_projective()) &&

           !C::is_on_curve(self.w_0.to_projective()) &&
           !C::is_on_curve(self.w_1.to_projective())
        
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
        writer.write_all(&self.bid_cm.x.to_bytes_le())?;
        writer.write_all(&self.bid_cm.y.to_bytes_le())?;

        writer.write_all(&self.r_cm.x.to_bytes_le())?;
        writer.write_all(&self.r_cm.y.to_bytes_le())?;

        writer.write_all(&self.r_inv_cm.x.to_bytes_le())?;
        writer.write_all(&self.r_inv_cm.y.to_bytes_le())?;

        writer.write_all(&self.f_cm.x.to_bytes_le())?;
        writer.write_all(&self.f_cm.y.to_bytes_le())?;

        writer.write_all(&self.diff_cm.x.to_bytes_le())?;
        writer.write_all(&self.diff_cm.y.to_bytes_le())?;

        writer.write_all(&self.g_cm.x.to_bytes_le())?;
        writer.write_all(&self.g_cm.y.to_bytes_le())?;

        writer.write_all(&self.q_price_opening.to_bytes_le())?;
        writer.write_all(&self.bid_opening.to_bytes_le())?;
        writer.write_all(&self.bid_shift_opening.to_bytes_le())?;
        writer.write_all(&self.f_opening.to_bytes_le())?;
        writer.write_all(&self.r_opening.to_bytes_le())?;
        writer.write_all(&self.r_inv_opening.to_bytes_le())?;
        writer.write_all(&self.diff_opening.to_bytes_le())?;
        writer.write_all(&self.g_opening.to_bytes_le())?;
        
        writer.write_all(&self.q_chunk_0_cm.x.to_bytes_le())?;
        writer.write_all(&self.q_chunk_0_cm.y.to_bytes_le())?;
        writer.write_all(&self.q_chunk_1_cm.x.to_bytes_le())?;
        writer.write_all(&self.q_chunk_1_cm.y.to_bytes_le())?;

        writer.write_all(&self.q_chunk_0_opening.to_bytes_le())?;
        writer.write_all(&self.q_chunk_1_opening.to_bytes_le())?;

        writer.write_all(&self.w_0.x.to_bytes_le())?;
        writer.write_all(&self.w_0.y.to_bytes_le())?;
        writer.write_all(&self.w_1.x.to_bytes_le())?;
        writer.write_all(&self.w_1.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        10 * scalar_len + 2 * base_len * 10
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

        let bid_cm = read_affine(&mut reader)?;
        let r_cm = read_affine(&mut reader)?;
        let r_inv_cm = read_affine(&mut reader)?;
        let f_cm = read_affine(&mut reader)?;
        let diff_cm = read_affine(&mut reader)?;
        let g_cm = read_affine(&mut reader)?;

        let mut q_price_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_price_opening_buf)?;
        let q_price_opening = C::ScalarField::from_bytes_le(&q_price_opening_buf);

        let mut bid_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut bid_opening_buf)?;
        let bid_opening = C::ScalarField::from_bytes_le(&bid_opening_buf);

        let mut bid_shift_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut bid_shift_opening_buf)?;
        let bid_shift_opening = C::ScalarField::from_bytes_le(&bid_shift_opening_buf);

        let mut f_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut f_opening_buf)?;
        let f_opening = C::ScalarField::from_bytes_le(&f_opening_buf);

        let mut r_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut r_opening_buf)?;
        let r_opening = C::ScalarField::from_bytes_le(&r_opening_buf);

        let mut r_inv_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut r_inv_opening_buf)?;
        let r_inv_opening = C::ScalarField::from_bytes_le(&r_inv_opening_buf);

        let mut diff_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut diff_opening_buf)?;
        let diff_opening = C::ScalarField::from_bytes_le(&diff_opening_buf);

        let mut g_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut g_opening_buf)?;
        let g_opening = C::ScalarField::from_bytes_le(&g_opening_buf);

        let q_chunk_0_cm = read_affine(&mut reader)?;
        let q_chunk_1_cm = read_affine(&mut reader)?;

        let mut q_chunk_0_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_chunk_0_opening_buf)?;
        let q_chunk_0_opening = C::ScalarField::from_bytes_le(&q_chunk_0_opening_buf);

        let mut q_chunk_1_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_chunk_1_opening_buf)?;
        let q_chunk_1_opening = C::ScalarField::from_bytes_le(&q_chunk_1_opening_buf);

        let w_0 = read_affine(&mut reader)?;
        let w_1 = read_affine(&mut reader)?;

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(bid_cm.to_projective()) &&
            !C::is_on_curve(r_cm.to_projective()) &&
            !C::is_on_curve(r_inv_cm.to_projective()) &&
            !C::is_on_curve(f_cm.to_projective()) &&
            !C::is_on_curve(diff_cm.to_projective()) &&
            !C::is_on_curve(g_cm.to_projective()) &&

            !C::is_on_curve(q_chunk_0_cm.to_projective()) &&
            !C::is_on_curve(q_chunk_1_cm.to_projective()) &&

            !C::is_on_curve(w_0.to_projective()) &&
            !C::is_on_curve(w_1.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { bid_cm, r_cm, r_inv_cm, f_cm, diff_cm, g_cm,
                  q_price_opening, bid_opening, bid_shift_opening, f_opening, r_opening, r_inv_opening, diff_opening, g_opening,
                  q_chunk_0_cm, q_chunk_1_cm,
                  q_chunk_0_opening, q_chunk_1_opening,
                  w_0, w_1
        })
    }
}