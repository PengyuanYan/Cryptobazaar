// This file contains the Structs for Gates and corresponding serlization functions.
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

// Check if data is valid
impl<C: Curve> Valid for VerifierIndex<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.q_price_cm.to_projective()) {
            return Err(SerializationError::InvalidData);
        }
        Ok(())
    }
}

// Serializition function for VerifierIndex
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

// Derializition function for VerifierIndex
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
    pub random_x: P,
    pub random_r: P,
    pub random_r_inv: P,
    pub diff_f: P,
    pub hidden_bid: P,
    pub e: PhantomData<C>,
}

pub struct Proof<C: Curve> {
    pub bid_cm: Affine::<C>,
    pub random_r_cm: Affine::<C>,
    pub random_r_inv_cm: Affine::<C>,
    pub random_x_cm: Affine::<C>,
    pub diff_f_cm: Affine::<C>,
    pub hidden_bid_cm: Affine::<C>,

    pub q_price_opening: C::ScalarField,
    pub bid_opening: C::ScalarField,
    pub bid_shift_opening: C::ScalarField,
    pub random_x_opening: C::ScalarField,
    pub random_r_opening: C::ScalarField,
    pub random_r_inv_opening: C::ScalarField,
    pub diff_f_opening: C::ScalarField,
    pub hidden_bid_opening: C::ScalarField,

    pub q_chunk_0_cm: Affine::<C>,
    pub q_chunk_1_cm: Affine::<C>,

    pub q_chunk_0_opening: C::ScalarField,
    pub q_chunk_1_opening: C::ScalarField,

    pub w_0: Affine::<C>,
    pub w_1: Affine::<C>,
}

// Check if data is valid
impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.bid_cm.to_projective()) &&
           !C::is_on_curve(self.random_r_cm.to_projective()) &&
           !C::is_on_curve(self.random_r_inv_cm.to_projective()) &&
           !C::is_on_curve(self.random_x_cm.to_projective()) &&
           !C::is_on_curve(self.diff_f_cm.to_projective()) &&
           !C::is_on_curve(self.hidden_bid_cm.to_projective()) &&

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

// Serializition function for Proof
impl<C: Curve> CanonicalSerialize for Proof<C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&self.bid_cm.x.to_bytes_le())?;
        writer.write_all(&self.bid_cm.y.to_bytes_le())?;

        writer.write_all(&self.random_r_cm.x.to_bytes_le())?;
        writer.write_all(&self.random_r_cm.y.to_bytes_le())?;

        writer.write_all(&self.random_r_inv_cm.x.to_bytes_le())?;
        writer.write_all(&self.random_r_inv_cm.y.to_bytes_le())?;

        writer.write_all(&self.random_x_cm.x.to_bytes_le())?;
        writer.write_all(&self.random_x_cm.y.to_bytes_le())?;

        writer.write_all(&self.diff_f_cm.x.to_bytes_le())?;
        writer.write_all(&self.diff_f_cm.y.to_bytes_le())?;

        writer.write_all(&self.hidden_bid_cm.x.to_bytes_le())?;
        writer.write_all(&self.hidden_bid_cm.y.to_bytes_le())?;

        writer.write_all(&self.q_price_opening.to_bytes_le())?;
        writer.write_all(&self.bid_opening.to_bytes_le())?;
        writer.write_all(&self.bid_shift_opening.to_bytes_le())?;
        writer.write_all(&self.random_x_opening.to_bytes_le())?;
        writer.write_all(&self.random_r_opening.to_bytes_le())?;
        writer.write_all(&self.random_r_inv_opening.to_bytes_le())?;
        writer.write_all(&self.diff_f_opening.to_bytes_le())?;
        writer.write_all(&self.hidden_bid_opening.to_bytes_le())?;
        
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

// Derializition function for Proof
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
        let random_r_cm = read_affine(&mut reader)?;
        let random_r_inv_cm = read_affine(&mut reader)?;
        let random_x_cm = read_affine(&mut reader)?;
        let diff_f_cm = read_affine(&mut reader)?;
        let hidden_bid_cm = read_affine(&mut reader)?;

        let mut q_price_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_price_opening_buf)?;
        let q_price_opening = C::ScalarField::from_bytes_le(&q_price_opening_buf);

        let mut bid_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut bid_opening_buf)?;
        let bid_opening = C::ScalarField::from_bytes_le(&bid_opening_buf);

        let mut bid_shift_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut bid_shift_opening_buf)?;
        let bid_shift_opening = C::ScalarField::from_bytes_le(&bid_shift_opening_buf);

        let mut random_x_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut random_x_opening_buf)?;
        let random_x_opening = C::ScalarField::from_bytes_le(&random_x_opening_buf);

        let mut random_r_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut random_r_opening_buf)?;
        let random_r_opening = C::ScalarField::from_bytes_le(&random_r_opening_buf);

        let mut random_r_inv_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut random_r_inv_opening_buf)?;
        let random_r_inv_opening = C::ScalarField::from_bytes_le(&random_r_inv_opening_buf);

        let mut diff_f_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut diff_f_opening_buf)?;
        let diff_f_opening = C::ScalarField::from_bytes_le(&diff_f_opening_buf);

        let mut hidden_bid_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut hidden_bid_opening_buf)?;
        let hidden_bid_opening = C::ScalarField::from_bytes_le(&hidden_bid_opening_buf);

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
            !C::is_on_curve(random_r_cm.to_projective()) &&
            !C::is_on_curve(random_r_inv_cm.to_projective()) &&
            !C::is_on_curve(random_x_cm.to_projective()) &&
            !C::is_on_curve(diff_f_cm.to_projective()) &&
            !C::is_on_curve(hidden_bid_cm.to_projective()) &&

            !C::is_on_curve(q_chunk_0_cm.to_projective()) &&
            !C::is_on_curve(q_chunk_1_cm.to_projective()) &&

            !C::is_on_curve(w_0.to_projective()) &&
            !C::is_on_curve(w_1.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { bid_cm, random_r_cm, random_r_inv_cm, random_x_cm, diff_f_cm, hidden_bid_cm,
                  q_price_opening, bid_opening, bid_shift_opening, random_x_opening, random_r_opening, random_r_inv_opening, diff_f_opening, hidden_bid_opening,
                  q_chunk_0_cm, q_chunk_1_cm,
                  q_chunk_0_opening, q_chunk_1_opening,
                  w_0, w_1
        })
    }
}

#[cfg(test)]
mod serialize_test {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_core::curve::Curve;
    use icicle_bn254::curve::ScalarCfg;
    
    use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, Compress};
    use ark_serialize::Validate;
    use icicle_core::traits::GenerateRandom;

    use super::{VerifierIndex, Proof};

    #[test]
    fn test_serialize() {
        let verifierindex = VerifierIndex::<Bn254CurveCfg> { 
                                q_price_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
        };
        
        let mut data = Vec::new();
        verifierindex.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = VerifierIndex::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(verifierindex.q_price_cm, result.q_price_cm);

        let proof = Proof::<Bn254CurveCfg> {
                            bid_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            random_r_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            random_r_inv_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            random_x_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            diff_f_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            hidden_bid_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],

                            q_price_opening: ScalarCfg::generate_random(1)[0],
                            bid_opening: ScalarCfg::generate_random(1)[0],
                            bid_shift_opening: ScalarCfg::generate_random(1)[0],
                            random_x_opening: ScalarCfg::generate_random(1)[0],
                            random_r_opening: ScalarCfg::generate_random(1)[0],
                            random_r_inv_opening: ScalarCfg::generate_random(1)[0],
                            diff_f_opening: ScalarCfg::generate_random(1)[0],
                            hidden_bid_opening: ScalarCfg::generate_random(1)[0],

                            q_chunk_0_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            q_chunk_1_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],

                            q_chunk_0_opening: ScalarCfg::generate_random(1)[0],
                            q_chunk_1_opening: ScalarCfg::generate_random(1)[0],
 
                            w_0: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            w_1: Bn254CurveCfg::generate_random_affine_points(1)[0],
        };

        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(proof.bid_cm, result.bid_cm);
        assert_eq!(proof.random_r_cm, result.random_r_cm);
        assert_eq!(proof.random_r_inv_cm, result.random_r_inv_cm);
        assert_eq!(proof.random_x_cm, result.random_x_cm);
        assert_eq!(proof.diff_f_cm, result.diff_f_cm);
        assert_eq!(proof.hidden_bid_cm, result.hidden_bid_cm);
        assert_eq!(proof.q_price_opening, result.q_price_opening);
        assert_eq!(proof.bid_opening, result.bid_opening);
        assert_eq!(proof.bid_shift_opening, result.bid_shift_opening);
        assert_eq!(proof.random_x_opening, result.random_x_opening);
        assert_eq!(proof.random_r_opening, result.random_r_opening);
        assert_eq!(proof.random_r_inv_opening, result.random_r_inv_opening);
        assert_eq!(proof.diff_f_opening, result.diff_f_opening);
        assert_eq!(proof.hidden_bid_opening, result.hidden_bid_opening);
        assert_eq!(proof.q_chunk_0_cm, result.q_chunk_0_cm);
        assert_eq!(proof.q_chunk_1_cm, result.q_chunk_1_cm);
        assert_eq!(proof.q_chunk_0_opening, result.q_chunk_0_opening);
        assert_eq!(proof.q_chunk_1_opening, result.q_chunk_1_opening);
        assert_eq!(proof.w_0, result.w_0);
        assert_eq!(proof.w_1, result.w_1);
    }
}