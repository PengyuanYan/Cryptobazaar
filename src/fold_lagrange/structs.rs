use icicle_core::traits::FieldImpl;
use icicle_core::curve::{Curve,Affine};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;
use ark_serialize::Valid;
use ark_serialize::SerializationError;

use crate::acc::structs::Proof as AccProof;
use crate::univariate_sumcheck::structs::Proof as UVProof;

use crate::{
    acc::structs::Error as AccError, univariate_sumcheck::structs::Error as SumcheckError,
};

#[derive(Debug)]
pub enum Error {
    AccError,
    SumcheckError,
}

impl From<AccError> for Error {
    fn from(_: AccError) -> Self {
        return Self::AccError;
    }
}

impl From<SumcheckError> for Error {
    fn from(_: SumcheckError) -> Self {
        return Self::SumcheckError;
    }
}

pub struct Instance<const N: usize, const LOG_N: usize, C: Curve> {
    pub lb_commitments: [Affine::<C>; N],
    pub challenges: [C::ScalarField; LOG_N],
}

pub struct Proof<C: Curve> {
    pub p_cm: Affine::<C>,
    pub acc_cm: Affine::<C>,
    pub acc_proof: AccProof<C>,
    pub sumcheck_proof: UVProof<C>,
}

impl<C: Curve> Valid for Proof<C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.p_cm.to_projective()) &&
           !C::is_on_curve(self.acc_cm.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        self.acc_proof.check()?;
        self.sumcheck_proof.check()?;

        Ok(())
    }
}

impl<C: Curve> CanonicalSerialize for Proof<C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&self.p_cm.x.to_bytes_le())?;
        writer.write_all(&self.p_cm.y.to_bytes_le())?;
        writer.write_all(&self.acc_cm.x.to_bytes_le())?;
        writer.write_all(&self.acc_cm.y.to_bytes_le())?;

        self.acc_proof.serialize_with_mode(&mut writer, _compress).unwrap();
        self.sumcheck_proof.serialize_with_mode(&mut writer, _compress).unwrap();

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();

        2 * base_len * 2 + self.acc_proof.serialized_size(_compress) + self.sumcheck_proof.serialized_size(_compress)
    }
}

impl<C: Curve> CanonicalDeserialize for Proof<C> {
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
        
        let p_cm = read_affine(&mut reader)?;
        let acc_cm = read_affine(&mut reader)?;
        
        let acc_proof = AccProof::<C>::deserialize_with_mode(&mut reader, _compress, validate).unwrap();
        let sumcheck_proof = UVProof::<C>::deserialize_with_mode(&mut reader, _compress, validate).unwrap();

        if matches!(validate, Validate::Yes) {
            if !C::is_on_curve(p_cm.to_projective()) ||
               !C::is_on_curve(acc_cm.to_projective())
            {
                return Err(SerializationError::InvalidData);
            }
        }

        Ok(Self { p_cm,
                  acc_cm,
                  acc_proof,
                  sumcheck_proof,
        })
    }
}

#[cfg(test)]
mod serialize_test {
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_core::curve::Curve;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;
    
    use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, Compress};
    use ark_serialize::Validate;

    use super::Proof;

    use crate::acc::structs::Proof as AccProof;
    use crate::univariate_sumcheck::structs::Proof as UVProof;

    #[test]
    fn test_serialize() {
        let acc_proof = AccProof::<Bn254CurveCfg> { 
                                  q: ScalarCfg::generate_random(1)[0],
                                  acc_opening: ScalarCfg::generate_random(1)[0],
                                  acc_shifted_opening: ScalarCfg::generate_random(1)[0],

                                  q_0: Bn254CurveCfg::generate_random_affine_points(1)[0],
                                  q_1: Bn254CurveCfg::generate_random_affine_points(1)[0],
                                  q_2: Bn254CurveCfg::generate_random_affine_points(1)[0],
        };
        
        let sumcheck_proof = UVProof::<Bn254CurveCfg> {
                                      r_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                                      q_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
    
                                      a_opening: ScalarCfg::generate_random(1)[0],
                                      b_opening: ScalarCfg::generate_random(1)[0],
                                      r_opening: ScalarCfg::generate_random(1)[0],
                                      q_opening: ScalarCfg::generate_random(1)[0],

                                      batch_opening_proof: Bn254CurveCfg::generate_random_affine_points(1)[0],
        };

        let proof = Proof::<Bn254CurveCfg> {
                         p_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                         acc_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                         acc_proof,
                         sumcheck_proof,
        };

        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.p_cm, proof.p_cm);
        assert_eq!(result.acc_cm, proof.acc_cm);
        assert_eq!(result.acc_proof.q, proof.acc_proof.q);
        assert_eq!(result.acc_proof.q_0, proof.acc_proof.q_0);
        assert_eq!(result.sumcheck_proof.q_cm, proof.sumcheck_proof.q_cm);
        assert_eq!(result.sumcheck_proof.batch_opening_proof, proof.sumcheck_proof.batch_opening_proof);
    }
}