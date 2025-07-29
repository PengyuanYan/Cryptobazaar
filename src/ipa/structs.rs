use icicle_core::traits::FieldImpl;
use icicle_core::curve::{Curve,Affine};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;
use ark_serialize::Valid;
use ark_serialize::SerializationError;

use crate::verifiable_folding_sumcheck::structs::Proof as VFSProof;

pub struct Instance<const N: usize, C: Curve> {
    pub ac: Affine::<C>,
    pub b: [Affine::<C>; N],
    pub h_base: Affine::<C>,
    pub c: [Affine::<C>; N],
}

impl<const N: usize, C: Curve> Valid for Instance<N, C> {
    fn check(&self) -> Result<(), SerializationError> {
        if !C::is_on_curve(self.ac.to_projective()) &&
           !C::is_on_curve(self.h_base.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }
        
        for i in self.b {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }

        for i in self.c {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }

        Ok(())
    }
}

impl<const N: usize, C: Curve> CanonicalSerialize for Instance<N, C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&self.ac.x.to_bytes_le())?;
        writer.write_all(&self.ac.y.to_bytes_le())?;
        
        writer.write_all(&(self.b.len() as u64).to_le_bytes())?;

        for i in self.b {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        writer.write_all(&self.h_base.x.to_bytes_le())?;
        writer.write_all(&self.h_base.y.to_bytes_le())?;

        for i in self.c {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        let vec_len = self.b.len();

        2 * base_len * (2 + 2 * vec_len) + 8
    }
}

impl<const N: usize, C: Curve> CanonicalDeserialize for Instance<N, C> {
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
        
        let ac = read_affine(&mut reader)?;
        
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let len_u64 = u64::from_le_bytes(len_buf);
        let len_usize = usize::try_from(len_u64)
            .map_err(|_| SerializationError::InvalidData)?;
        
        let mut b = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let b_i = read_affine(&mut reader)?;
            b.push(b_i);
        }

        let h_base = read_affine(&mut reader)?;

        let mut c = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let c_i = read_affine(&mut reader)?;
            c.push(c_i);
        }

        if matches!(validate, Validate::Yes) &&
            !C::is_on_curve(ac.to_projective()) &&
            !C::is_on_curve(h_base.to_projective())
        {
            return Err(SerializationError::InvalidData);
        }

        Ok(Self { ac,
                  b: b.try_into().unwrap(),
                  h_base,
                  c: c.try_into().unwrap(),
        })
    }
}

pub struct Witness<const N: usize, C: Curve> {
    pub a: [C::ScalarField; N],
}

pub struct Proof<const LOG_N: usize, C: Curve> {
    pub l: [Affine::<C>; LOG_N],
    pub r: [Affine::<C>; LOG_N],

    pub vfs_proof: VFSProof<C>,
}

impl<const LOG_N: usize, C: Curve> Valid for Proof<LOG_N, C> {
    fn check(&self) -> Result<(), SerializationError> {
        
        self.vfs_proof.check()?;

        for i in self.l {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }

        for i in self.r {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }


        Ok(())
    }
}

impl<const LOG_N: usize, C: Curve> CanonicalSerialize for Proof<LOG_N, C> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress, // it cant change anything
    ) -> Result<(), SerializationError> {
        writer.write_all(&(self.l.len() as u64).to_le_bytes())?;

        for i in self.l {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        for i in self.r {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }
        
        self.vfs_proof.serialize_with_mode(&mut writer, _compress).unwrap();

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        let vec_len = self.l.len();

        2 * base_len * (2 * vec_len) + 8 + self.vfs_proof.serialized_size(_compress)
    }
}

impl<const LOG_N: usize, C: Curve> CanonicalDeserialize for Proof<LOG_N, C> {
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
        
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let len_u64 = u64::from_le_bytes(len_buf);
        let len_usize = usize::try_from(len_u64)
            .map_err(|_| SerializationError::InvalidData)?;
        
        let mut l = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let l_i = read_affine(&mut reader)?;
            l.push(l_i);
        }

        let mut r = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let r_i = read_affine(&mut reader)?;
            r.push(r_i);
        }
        
        let vfs_proof = VFSProof::<C>::deserialize_with_mode(reader, _compress, validate).unwrap();

        if matches!(validate, Validate::Yes) {
            if l.iter().any(|pt| !C::is_on_curve(pt.to_projective())) ||
               r.iter().any(|pt| !C::is_on_curve(pt.to_projective()))
            {
                return Err(SerializationError::InvalidData);
            }
        }

        Ok(Self { l: l.try_into().unwrap(),
                  r: r.try_into().unwrap(),
                  vfs_proof,
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

    use std::array::from_fn;

    use crate::verifiable_folding_sumcheck::structs::Proof as VFSProof;
    
    const N: usize = 4;
    const LOG_N: usize = 2;

    #[test]
    fn test_serialize() {
        let instance = Instance::<N, Bn254CurveCfg> { 
                                ac: Affine::<Bn254CurveCfg>::zero(),
                                b: from_fn(|_| Affine::<Bn254CurveCfg>::zero()),
                                h_base: Affine::<Bn254CurveCfg>::zero(),
                                c: from_fn(|_| Affine::<Bn254CurveCfg>::zero()),
        };
        
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Instance::<N, Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(instance.b.len(), result.b.len());
        assert_eq!(instance.c.len(), result.c.len());
        assert_eq!(instance.c.len(), N);

        let vfs_proof = VFSProof::<Bn254CurveCfg> { 
                                  s: Affine::<Bn254CurveCfg>::zero(),
                                  blinder_cm: Affine::<Bn254CurveCfg>::zero(),

                                  z_1: Bn254ScalarField::one(),
                                  z_2: Bn254ScalarField::one() + Bn254ScalarField::one(),
                                  r_cm: Affine::<Bn254CurveCfg>::zero(),
                                  r_degree_cm: Affine::<Bn254CurveCfg>::zero(),
                                  q_cm: Affine::<Bn254CurveCfg>::zero(),

                                  a_opening: Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one(),
                                  blinder_opening: Bn254ScalarField::one(),
                                  r_opening: Bn254ScalarField::one() + Bn254ScalarField::one(),
                                  q_opening: Bn254ScalarField::one() + Bn254ScalarField::one() + Bn254ScalarField::one(),

                                  batch_opening_proof: Affine::<Bn254CurveCfg>::zero(),
                                };
        
        let proof = Proof::<LOG_N, Bn254CurveCfg> {
                            l: from_fn(|_| Affine::<Bn254CurveCfg>::zero()),
                            r: from_fn(|_| Affine::<Bn254CurveCfg>::zero()),
                            vfs_proof
        };

        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<LOG_N, Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(result.l.len(), proof.l.len());
        assert_eq!(result.r.len(), proof.r.len());
        assert_eq!(result.vfs_proof.z_2, proof.vfs_proof.z_2);
        assert_eq!(result.vfs_proof.r_opening, proof.vfs_proof.r_opening);
    }
}