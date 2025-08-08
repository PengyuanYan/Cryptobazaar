use icicle_core::traits::FieldImpl;
use icicle_core::curve::{Curve,Affine};
use icicle_core::polynomials::UnivariatePolynomial;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;
use ark_serialize::Valid;
use ark_serialize::SerializationError;

#[derive(Debug)]
pub enum Error {
    OpeningFailed,
    RelationCheckFailed,
    PedersenOpeningFailed,
    DegreeCheckShiftMissing,
    DegreeCheckFailed,
}

pub struct Instance<C: Curve> {
    pub n: usize,
    pub p_base: Affine::<C>,
    pub h_base: Affine::<C>,
    pub a_cm: Affine::<C>,
    pub pedersen: Affine::<C>,
    pub challenges: Vec<C::ScalarField>,
}

impl<C: Curve> Valid for Instance<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.p_base.to_projective()) &&
           !C::is_on_curve(self.h_base.to_projective()) &&
           !C::is_on_curve(self.a_cm.to_projective()) &&
           !C::is_on_curve(self.pedersen.to_projective())
        
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

        writer.write_all(&self.p_base.x.to_bytes_le())?;
        writer.write_all(&self.p_base.y.to_bytes_le())?;

        writer.write_all(&self.h_base.x.to_bytes_le())?;
        writer.write_all(&self.h_base.y.to_bytes_le())?;

        writer.write_all(&self.a_cm.x.to_bytes_le())?;
        writer.write_all(&self.a_cm.y.to_bytes_le())?;

        writer.write_all(&self.pedersen.x.to_bytes_le())?;
        writer.write_all(&self.pedersen.y.to_bytes_le())?;

        writer.write_all(&(self.challenges.len() as u64).to_le_bytes())?;

        for i in 0..self.challenges.len() {
            writer.write_all(&self.challenges[i].to_bytes_le())?;
        }

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let len = self.challenges.len();
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        16 + len * scalar_len + 2 * base_len * 4
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
        
        let p_base = read_affine(&mut reader)?;
        let h_base = read_affine(&mut reader)?;
        let a_cm = read_affine(&mut reader)?;
        let pedersen = read_affine(&mut reader)?;

        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let len_buf_u64 = u64::from_le_bytes(len_buf);
        let len_usize = usize::try_from(len_buf_u64)
            .map_err(|_| SerializationError::InvalidData)?;

        let mut challenges = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let mut challenge_buf = vec![0u8; scalar_len];
            reader.read_exact(&mut challenge_buf)?;
            let challenge = C::ScalarField::from_bytes_le(&challenge_buf);
            challenges.push(challenge);
        }

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(p_base.to_projective()) &&
            !C::is_on_curve(h_base.to_projective()) &&
            !C::is_on_curve(a_cm.to_projective()) &&
            !C::is_on_curve(pedersen.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { n: n_usize,
                  p_base, h_base, a_cm, pedersen,
                  challenges,
        })
    }
}

pub struct Witness<C, P> 
where
    C: Curve,
    P: UnivariatePolynomial<Field = C::ScalarField>,
{
    pub a: P,
    pub x: C::ScalarField,
    pub r: C::ScalarField,
}

pub struct Proof<C: Curve> {
    // round 1
    pub(crate) s: Affine::<C>,
    pub(crate) blinder_cm: Affine::<C>,

    // round 2
    pub(crate) z_1: C::ScalarField,
    pub(crate) z_2: C::ScalarField,
    pub(crate) r_cm: Affine::<C>,
    pub(crate) r_degree_cm: Affine::<C>,
    pub(crate) q_cm: Affine::<C>,

    // round 3
    pub(crate) a_opening: C::ScalarField,
    pub(crate) blinder_opening: C::ScalarField,
    pub(crate) r_opening: C::ScalarField,
    pub(crate) q_opening: C::ScalarField,

    // round 4
    pub(crate) batch_opening_proof: Affine::<C>,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.s.to_projective()) &&
           !C::is_on_curve(self.blinder_cm.to_projective()) &&
           !C::is_on_curve(self.r_cm.to_projective()) &&
           !C::is_on_curve(self.r_degree_cm.to_projective()) &&
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
        writer.write_all(&self.s.x.to_bytes_le())?;
        writer.write_all(&self.s.y.to_bytes_le())?;
        writer.write_all(&self.blinder_cm.x.to_bytes_le())?;
        writer.write_all(&self.blinder_cm.y.to_bytes_le())?;

        writer.write_all(&self.z_1.to_bytes_le())?;
        writer.write_all(&self.z_2.to_bytes_le())?;

        writer.write_all(&self.r_cm.x.to_bytes_le())?;
        writer.write_all(&self.r_cm.y.to_bytes_le())?;
        writer.write_all(&self.r_degree_cm.x.to_bytes_le())?;
        writer.write_all(&self.r_degree_cm.y.to_bytes_le())?;
        writer.write_all(&self.q_cm.x.to_bytes_le())?;
        writer.write_all(&self.q_cm.y.to_bytes_le())?;

        writer.write_all(&self.a_opening.to_bytes_le())?;
        writer.write_all(&self.blinder_opening.to_bytes_le())?;
        writer.write_all(&self.r_opening.to_bytes_le())?;
        writer.write_all(&self.q_opening.to_bytes_le())?;

        writer.write_all(&self.batch_opening_proof.x.to_bytes_le())?;
        writer.write_all(&self.batch_opening_proof.y.to_bytes_le())?;

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let scalar_len = C::ScalarField::zero().to_bytes_le().len();
        let base_len = C::BaseField::zero().to_bytes_le().len();

        6 * scalar_len + 2 * base_len * 6
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
        
        let s = read_affine(&mut reader)?;
        let blinder_cm = read_affine(&mut reader)?;

        let mut z_1_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_1_buf)?;
        let z_1 = C::ScalarField::from_bytes_le(&z_1_buf);
        let mut z_2_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut z_2_buf)?;
        let z_2 = C::ScalarField::from_bytes_le(&z_2_buf);

        let r_cm = read_affine(&mut reader)?;
        let r_degree_cm = read_affine(&mut reader)?;
        let q_cm = read_affine(&mut reader)?;

        let mut a_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut a_opening_buf)?;
        let a_opening = C::ScalarField::from_bytes_le(&a_opening_buf);
        let mut blinder_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut blinder_opening_buf)?;
        let blinder_opening = C::ScalarField::from_bytes_le(&blinder_opening_buf);
        let mut r_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut r_opening_buf)?;
        let r_opening = C::ScalarField::from_bytes_le(&r_opening_buf);
        let mut q_opening_buf = vec![0u8; scalar_len];
        reader.read_exact(&mut q_opening_buf)?;
        let q_opening = C::ScalarField::from_bytes_le(&q_opening_buf);

        let batch_opening_proof = read_affine(&mut reader)?;

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(s.to_projective()) &&
            !C::is_on_curve(blinder_cm.to_projective()) &&
            !C::is_on_curve(r_cm.to_projective()) &&
            !C::is_on_curve(r_degree_cm.to_projective()) &&
            !C::is_on_curve(q_cm.to_projective()) &&
            !C::is_on_curve(batch_opening_proof.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { s, blinder_cm,
                  z_1, z_2, r_cm, r_degree_cm, q_cm,
                  a_opening, blinder_opening, r_opening, q_opening,
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
        let p_base = Affine::<Bn254CurveCfg>::zero();
        let h_base = Affine::<Bn254CurveCfg>::zero();
        let a_cm = Affine::<Bn254CurveCfg>::zero();
        let pedersen = Affine::<Bn254CurveCfg>::zero();
        let challenges = vec![(Bn254ScalarField::one() + Bn254ScalarField::one()), (Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one())];

        let instance = Instance { n,
                                  p_base, h_base, a_cm, pedersen,
                                  challenges,
                                };
        
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Instance::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(instance.n, result.n);
        assert_eq!(instance.challenges.len(), result.challenges.len());
        assert_eq!(instance.challenges[0], result.challenges[0]);
        assert_eq!(instance.challenges[1], result.challenges[1]);

        let s = Affine::<Bn254CurveCfg>::zero();
        let blinder_cm = Affine::<Bn254CurveCfg>::zero();
        
        let z_1 = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let z_2 = Bn254ScalarField::one();
        let r_cm = Affine::<Bn254CurveCfg>::zero();
        let r_degree_cm = Affine::<Bn254CurveCfg>::zero();
        let q_cm = Affine::<Bn254CurveCfg>::zero();

        let a_opening = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();
        let blinder_opening = Bn254ScalarField::one();
        let r_opening = Bn254ScalarField::one() + Bn254ScalarField::one();
        let q_opening = Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one();

        let batch_opening_proof = Affine::<Bn254CurveCfg>::zero();

        let proof = Proof { s, blinder_cm,
                            z_1, z_2, r_cm, r_degree_cm, q_cm,
                            a_opening, blinder_opening, r_opening, q_opening,
                            batch_opening_proof,
                           };
        
        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.z_1, proof.z_1);
        assert_eq!(result.blinder_opening, proof.blinder_opening);
    }
}