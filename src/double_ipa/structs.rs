use icicle_core::traits::FieldImpl;
use icicle_core::curve::{Curve,Affine};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::io::{Read, Write};
use ark_serialize::Validate;
use ark_serialize::Compress;
use ark_serialize::Valid;
use ark_serialize::SerializationError;

use crate::double_pedersen_schnorr::structs::Proof as PSProof;
use crate::fold_lagrange::structs::Proof as LFProof;

pub struct Instance<const N: usize, C: Curve> {
    pub ac: Affine::<C>,
    pub b: [Affine::<C>; N],
    pub h_base: Affine::<C>,
    pub c: [Affine::<C>; N],
    pub lagrange_basis: [Affine::<C>; N],
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

        for i in self.lagrange_basis {
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

        for i in self.lagrange_basis {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        let vec_len = self.b.len();

        2 * base_len * (2 + 3 * vec_len) + 8
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
        
        let mut lagrange_basis = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let lagrange_basis_i = read_affine(&mut reader)?;
            lagrange_basis.push(lagrange_basis_i);
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
                  lagrange_basis: lagrange_basis.try_into().unwrap(),
        })
    }
}

pub struct Witness<const N: usize, C: Curve> {
    pub a: [C::ScalarField; N],
}

pub struct Proof<const LOG_N: usize, C: Curve> {
    pub l_1: [Affine::<C>; LOG_N],
    pub r_1: [Affine::<C>; LOG_N],

    pub lf_proof: LFProof<C>,

    pub l_2: [Affine::<C>; LOG_N],
    pub r_2: [Affine::<C>; LOG_N],

    pub ps_proof: PSProof<C>,
}

impl<const LOG_N: usize, C: Curve> Valid for Proof<LOG_N, C> {
    fn check(&self) -> Result<(), SerializationError> {
        
        self.lf_proof.check()?;
        self.ps_proof.check()?;

        for i in self.l_1 {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }

        for i in self.r_1 {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }

        for i in self.l_2 {
            if !C::is_on_curve(i.to_projective()) {
                return Err(SerializationError::InvalidData);
            }
        }

        for i in self.l_2 {
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
        writer.write_all(&(self.l_1.len() as u64).to_le_bytes())?;

        for i in self.l_1 {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        for i in self.r_1 {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        for i in self.l_2 {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }

        for i in self.r_2 {
            writer.write_all(&i.x.to_bytes_le())?;
            writer.write_all(&i.y.to_bytes_le())?;
        }
        
        self.lf_proof.serialize_with_mode(&mut writer, _compress).unwrap();
        self.ps_proof.serialize_with_mode(&mut writer, _compress).unwrap();

        Ok(())
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        let base_len = C::BaseField::zero().to_bytes_le().len();
        let vec_len = self.l_1.len();

        2 * base_len * (2 * 2 * vec_len) + 8 + self.lf_proof.serialized_size(_compress) + self.ps_proof.serialized_size(_compress)
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
        
        let mut l_1 = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let l_i = read_affine(&mut reader)?;
            l_1.push(l_i);
        }

        let mut r_1 = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let r_i = read_affine(&mut reader)?;
            r_1.push(r_i);
        }

        let mut l_2 = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let l_i = read_affine(&mut reader)?;
            l_2.push(l_i);
        }

        let mut r_2 = Vec::with_capacity(len_usize);

        for _ in 0..len_usize {
            let r_i = read_affine(&mut reader)?;
            r_2.push(r_i);
        }
        
        let lf_proof = LFProof::<C>::deserialize_with_mode(&mut reader, _compress, validate).unwrap();
        let ps_proof = PSProof::<C>::deserialize_with_mode(&mut reader, _compress, validate).unwrap();

        if matches!(validate, Validate::Yes) {
            if l_1.iter().any(|pt| !C::is_on_curve(pt.to_projective())) ||
               r_1.iter().any(|pt| !C::is_on_curve(pt.to_projective())) ||
               l_2.iter().any(|pt| !C::is_on_curve(pt.to_projective())) ||
               r_2.iter().any(|pt| !C::is_on_curve(pt.to_projective()))
            {
                return Err(SerializationError::InvalidData);
            }
        }

        Ok(Self { l_1: l_1.try_into().unwrap(),
                  r_1: r_1.try_into().unwrap(),
                  lf_proof,
                  l_2: l_2.try_into().unwrap(),
                  r_2: r_2.try_into().unwrap(),
                  ps_proof
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

    use super::{Instance, Proof};

    use std::array::from_fn;

    use crate::double_pedersen_schnorr::structs::Proof as PSProof;
    use crate::fold_lagrange::structs::Proof as LFProof;
    use crate::acc::structs::Proof as AccProof;
    use crate::univariate_sumcheck::structs::Proof as UVProof;
    
    const N: usize = 4;
    const LOG_N: usize = 2;

    #[test]
    fn test_serialize() {
        let instance = Instance::<N, Bn254CurveCfg> { 
                                ac: Bn254CurveCfg::generate_random_affine_points(1)[0],
                                b: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
                                h_base: Bn254CurveCfg::generate_random_affine_points(1)[0],
                                c: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
                                lagrange_basis: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
        };
        
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Instance::<N, Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();
        
        assert_eq!(instance.b.len(), result.b.len());
        assert_eq!(instance.b[1], result.b[1]);
        assert_eq!(instance.c.len(), result.c.len());
        assert_eq!(instance.c[1], result.c[1]);
        assert_eq!(instance.c.len(), N);
        assert_eq!(instance.lagrange_basis.len(), result.lagrange_basis.len());
        assert_eq!(instance.lagrange_basis[1], result.lagrange_basis[1]);

        
        let psproof = PSProof::<Bn254CurveCfg> {
                               rand_1: Bn254CurveCfg::generate_random_affine_points(1)[0],
                               rand_2: Bn254CurveCfg::generate_random_affine_points(1)[0],
                            
                               z_1: ScalarCfg::generate_random(1)[0],
                               z_2: ScalarCfg::generate_random(1)[0],
                               z_3: ScalarCfg::generate_random(1)[0],

        };
        
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

        let lfproof = LFProof::<Bn254CurveCfg> {
                               p_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                               acc_cm: Bn254CurveCfg::generate_random_affine_points(1)[0],
                               acc_proof,
                               sumcheck_proof,
        };

        let proof = Proof::<LOG_N, Bn254CurveCfg> {
                           l_1: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
                           r_1: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
                           lf_proof: lfproof,

                           l_2: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
                           r_2: from_fn(|_| Bn254CurveCfg::generate_random_affine_points(1)[0]),
                           ps_proof: psproof
        };

        let mut data = Vec::new();
        proof.serialize_with_mode(&mut data, Compress::No).unwrap();

        let mut reader: &[u8] = &data;
        let result = Proof::<LOG_N, Bn254CurveCfg>::deserialize_with_mode(&mut reader, Compress::No, Validate::No).unwrap();

        assert_eq!(proof.l_1.len(), result.l_1.len());
        assert_eq!(proof.l_2.len(), result.l_2.len());
        assert_eq!(proof.r_1.len(), result.r_1.len());
        assert_eq!(proof.r_2.len(), result.r_2.len());
        assert_eq!(proof.l_1[0], result.l_1[0]);
        assert_eq!(proof.l_2[0], result.l_2[0]);
        assert_eq!(proof.r_1[0], result.r_1[0]);
        assert_eq!(proof.r_2[0], result.r_2[0]);

        assert_eq!(proof.lf_proof.p_cm, result.lf_proof.p_cm);
        assert_eq!(proof.lf_proof.acc_cm, result.lf_proof.acc_cm);
        assert_eq!(proof.lf_proof.acc_proof.acc_opening, result.lf_proof.acc_proof.acc_opening);
        assert_eq!(proof.lf_proof.acc_proof.q_1, result.lf_proof.acc_proof.q_1);
        assert_eq!(proof.lf_proof.sumcheck_proof.q_cm, result.lf_proof.sumcheck_proof.q_cm);
        assert_eq!(proof.lf_proof.sumcheck_proof.b_opening, result.lf_proof.sumcheck_proof.b_opening);
        assert_eq!(proof.ps_proof.rand_1, result.ps_proof.rand_1);
        assert_eq!(proof.ps_proof.z_3, result.ps_proof.z_3);
    }
}