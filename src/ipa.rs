use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::traits::FieldImpl;
use icicle_core::pairing::Pairing;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_runtime::memory::HostSlice;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, NTT, release_domain};
use icicle_core::traits::GenerateRandom;
use icicle_core::traits::Arithmetic;

use icicle_core::{msm, msm::MSMConfig};

use crate::kzg::{DegreeCheckVK, PK as KzgPk, VK as KzgVk};
use crate::utils::{msm_gpu, ntt_gpu};
use crate::utils::folding::{compute_folding_coeffs, AffFold, FFold, Fold};
use crate::verifiable_folding_sumcheck::{
    structs::{Error as VFSError, Instance as VFSInstance, Witness as VFSWitness},
    Argument as VFSArgument,
};

use self::structs::{Instance, Proof, Witness};
use self::tr::Transcript;

pub mod structs;
pub mod tr;

#[derive(Debug)]
pub enum Error {
    FoldShapeMismatch { round: usize, a_left: usize, b_right: usize, a_right: usize, b_left: usize },
    FoldWrongLength { a_len: usize, b_len: usize },
    SanityMismatch { position: usize },
}

pub struct InnerProduct<const N: usize, const LOG_N: usize, C1, C2, F, U>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    _e: PhantomData<(C1, C2, F, U)>,
}

impl<const N: usize, const LOG_N: usize, C1, C2, F, U> InnerProduct<N, LOG_N, C1, C2, F, U>
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    pub fn prove (
        instance: &Instance<N, C1>,
        witness: &Witness<N, C1>,
        pk: &KzgPk<C1, C2, F>,
        cpu_or_gpu: usize,
    ) -> Result<Proof<LOG_N, C1>, Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: GenerateRandom<<C1 as Curve>::ScalarField>,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField> + NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <U as UnivariatePolynomial>::FieldConfig: GenerateRandom<<C1 as Curve>::ScalarField>,
    {
        let mut acc_blinders = C1::ScalarField::zero();

        let mut l_msgs = Vec::<Affine::<C1>>::with_capacity(LOG_N);
        let mut r_msgs = Vec::<Affine::<C1>>::with_capacity(LOG_N);

        let mut alphas = Vec::with_capacity(LOG_N);

        let mut tr = Transcript::<N, LOG_N, C1>::new(b"ipa");
        tr.send_instance(instance);

        let r = tr.get_r();
        let r_pows: Vec<_> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * r)).take(N).collect();
        
        let cfg = MSMConfig::default();
        let mut c_fold_projective = if cpu_or_gpu == 0 {
            let mut c_fold_projective_v = vec![Projective::<C1>::zero(); 1];
            msm::msm(
                HostSlice::from_slice(&r_pows),
                HostSlice::from_slice(&instance.c),
                &cfg,
                HostSlice::from_mut_slice(&mut c_fold_projective_v),
            )
            .unwrap();

            c_fold_projective_v[0]
        } else {
             msm_gpu(&r_pows, &instance.c)
        };

        let mut a_folded = witness.a.clone().to_vec();

        let mut b_folded: Vec<Affine::<C1>> = Vec::with_capacity(instance.b.len());

        for i in 0..instance.b.len() {
            let bi = instance.b[i];
            let ri = r_pows[i];

            let result_projective = C1::mul_scalar(bi.to_projective(), ri);

            let mut result_affine = Affine::<C1>::zero();
            C1::to_affine(&result_projective, &mut result_affine);

            b_folded.push(result_affine);
        }
        
        for round in 0..(LOG_N) as usize {
            let a_left = &a_folded[..a_folded.len() / 2];
            let a_right = &a_folded[a_folded.len() / 2..];

            let b_left = &b_folded[..b_folded.len() / 2];
            let b_right = &b_folded[b_folded.len() / 2..];

            if a_left.len() != b_right.len() || a_right.len() != b_left.len() {
                return Err(Error::FoldShapeMismatch {
                    round: round,
                    a_left: a_left.len(),
                    b_right: b_right.len(),
                    a_right: a_right.len(),
                    b_left: b_left.len(),
                });
            }

            let l = if a_left.len() == 1 { 
                C1::mul_scalar(b_right[0].to_projective(), a_left[0])
            } else {
                let result = if cpu_or_gpu == 0 {
                    let mut buf = [Projective::<C1>::zero()];
                    msm::msm(
                        HostSlice::from_slice(a_left),
                        HostSlice::from_slice(b_right),
                        &cfg,
                        HostSlice::from_mut_slice(&mut buf),
                    )
                    .unwrap();

                    buf[0]
                } else {
                    msm_gpu(a_left, b_right)
                };

                result
            };

            let r = if a_right.len() == 1 { 
                C1::mul_scalar(b_left[0].to_projective(), a_right[0])
            } else {
                let result = if cpu_or_gpu == 0 {
                    let mut buf = [Projective::<C1>::zero()];
                    msm::msm(
                        HostSlice::from_slice(a_right),
                        HostSlice::from_slice(b_left),
                        &cfg,
                        HostSlice::from_mut_slice(&mut buf),
                    )
                    .unwrap();

                    buf[0]
                } else {
                    msm_gpu(a_right, b_left)
                };

                result
            };
            
            let blinders = <<C1::ScalarField as FieldImpl>::Config as GenerateRandom<C1::ScalarField>>::generate_random(2);
            let blinder_l = blinders[0];
            let blinder_r = blinders[1];

            let l_projective = l + C1::mul_scalar(instance.h_base.to_projective(), blinder_l);
            let mut l_affine = Affine::<C1>::zero();
            C1::to_affine(&l_projective, &mut l_affine);
            l_msgs.push(l_affine);

            let r_projective = r + C1::mul_scalar(instance.h_base.to_projective(), blinder_r);
            let mut r_affine = Affine::<C1>::zero();
            C1::to_affine(&r_projective, &mut r_affine);
            r_msgs.push(r_affine);

            tr.send_l_r(&l_affine, &r_affine);
            let alpha = tr.get_alpha_i();
            let alpha_inv = alpha.inv();
            alphas.push(alpha);

            acc_blinders = acc_blinders + alpha_inv * blinder_l + alpha * blinder_r;

            // fold vectors
            a_folded = FFold::<C1>::fold_vec(&a_folded, alpha).unwrap();
            b_folded = AffFold::<C1>::fold_vec(&b_folded, alpha_inv).unwrap();

            // derive new cm
            c_fold_projective = C1::mul_scalar(l_projective, alpha_inv) + c_fold_projective + C1::mul_scalar(r_projective, alpha);
        }
        
        // sanity
        if a_folded.len() != 1 || b_folded.len() != 1 {
            return Err(Error::FoldWrongLength {
                a_len: a_folded.len(),
                b_len: b_folded.len(),
            });
        }

        let lhs = C1::mul_scalar(b_folded[0].to_projective(), a_folded[0])
                + C1::mul_scalar(instance.h_base.to_projective(), acc_blinders);

        if lhs != c_fold_projective {
            return Err(Error::SanityMismatch {
                position: 1usize
            });
        }
        
        let mut c_fold_affine = Affine::<C1>::zero();
        C1::to_affine(&c_fold_projective, &mut c_fold_affine);
        
        let vfs_instance = VFSInstance::<C1> {
            n: N,
            p_base: b_folded[0],
            h_base: instance.h_base,
            a_cm: instance.ac,
            pedersen: c_fold_affine,
            challenges: alphas.try_into().unwrap(),
        };
        
        //let domain = get_root_of_unity::<C1::ScalarField>((N).try_into().unwrap());
        let cfg = NTTConfig::<C1::ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut a_coeffs = if cpu_or_gpu == 0 {
            let mut ntt_result =vec![C1::ScalarField::zero(); witness.a.len()];
            ntt(
                HostSlice::from_slice(&witness.a),
                NTTDir::kInverse,
                &cfg,
                HostSlice::from_mut_slice(&mut ntt_result),
            )
            .unwrap();

            ntt_result
        } else {
            let coset_gen = C1::ScalarField::one();
            ntt_gpu::<C1>(&witness.a, &coset_gen, NTTDir::kInverse)
        };
        
        //release_domain::<C1::ScalarField>().unwrap();

        let a = U::from_coeffs(HostSlice::from_slice(&a_coeffs), a_coeffs.len());

        let vfs_witness = VFSWitness {
            a,
            x: a_folded[0],
            r: acc_blinders,
        };
        
        // sanity
        let x = C1::mul_scalar(b_folded[0].to_projective(), a_folded[0])
              + C1::mul_scalar(instance.h_base.to_projective(), acc_blinders);

        if x != c_fold_projective {
            return Err(Error::SanityMismatch {
                position: 2usize
            });
        }

        let vfs_proof = VFSArgument::<C1, C2, F, U>::prove(&vfs_instance, &vfs_witness, pk);
        
        Ok(Proof {
            l: l_msgs.try_into().unwrap(),
            r: r_msgs.try_into().unwrap(),
            vfs_proof,
        })
    }

    pub fn verify(
        instance: &Instance<N, C1>,
        proof: &Proof<LOG_N, C1>,
        vk: &KzgVk<C1, C2, F>,
        degree_check_vk: &DegreeCheckVK<C1, C2, F>,
        cpu_or_gpu: usize,
    ) -> Result<(), VFSError>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <U as UnivariatePolynomial>::FieldConfig: GenerateRandom<<C1 as Curve>::ScalarField>,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        let mut tr = Transcript::<N, LOG_N, _>::new(b"ipa");
        tr.send_instance(instance);
        
        let r = tr.get_r();
        let r_pows: Vec<_> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * r)).take(N).collect();

        let cfg = MSMConfig::default();
        let mut c_fold_projective = if cpu_or_gpu == 0 {
            let mut c_fold_projective_v = vec![Projective::<C1>::zero(); 1];
            msm::msm(
                HostSlice::from_slice(&r_pows),
                HostSlice::from_slice(&instance.c),
                &cfg,
                HostSlice::from_mut_slice(&mut c_fold_projective_v),
            )
            .unwrap();

            c_fold_projective_v[0]
        } else {
             msm_gpu(&r_pows, &instance.c)
        };

        let mut b_rescaled: Vec<Affine::<C1>> = Vec::with_capacity(instance.b.len());

        for i in 0..instance.b.len() {
            let bi = instance.b[i];
            let ri = r_pows[i];

            let result_projective = C1::mul_scalar(bi.to_projective(), ri);

            let mut result_affine = Affine::<C1>::zero();
            C1::to_affine(&result_projective, &mut result_affine);

             b_rescaled.push(result_affine);
        }

        let mut alphas = Vec::with_capacity(LOG_N);
        let mut alpha_invs = Vec::with_capacity(LOG_N);

        for i in 0..LOG_N as usize {
            tr.send_l_r(&proof.l[i], &proof.r[i]);

            let alpha = tr.get_alpha_i();
            let alpha_inv = alpha.inv();
            alphas.push(alpha);
            alpha_invs.push(alpha_inv);

            c_fold_projective = C1::mul_scalar(proof.l[i].to_projective(), alpha_inv) + c_fold_projective + C1::mul_scalar(proof.r[i].to_projective(), alpha);
        }

        let mut c_fold_affine = Affine::<C1>::zero();
        C1::to_affine(&c_fold_projective, &mut c_fold_affine);

        // now fold the basis b and check pedersen openings
        let fold_coeffs = compute_folding_coeffs::<C1>(&alpha_invs);

        let mut b_folded_projective_v = vec![Projective::<C1>::zero(); 1];
        msm::msm(
            HostSlice::from_slice(&fold_coeffs),
            HostSlice::from_slice(&b_rescaled),
            &cfg,
            HostSlice::from_mut_slice(&mut b_folded_projective_v),
        )
        .unwrap();

        let b_folded_projective = b_folded_projective_v[0];
        let mut b_folded_affine = Affine::<C1>::zero();
        C1::to_affine(&b_folded_projective, &mut b_folded_affine);

        let vfs_instance = VFSInstance::<C1> {
            n: N,
            p_base: b_folded_affine,
            h_base: instance.h_base,
            a_cm: instance.ac,
            pedersen: c_fold_affine,
            challenges: alphas.try_into().unwrap(),
        };

        VFSArgument::<C1, C2, F, U>::verify(&vfs_instance, &proof.vfs_proof, vk, degree_check_vk)?;
        Ok(())
    }
}

#[cfg(test)]
mod ipa_tests {
    use std::collections::BTreeMap;

    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::traits::FieldImpl;
    use std::marker::PhantomData;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    use icicle_core::ntt::get_root_of_unity;
    use icicle_runtime::memory::HostSlice;
    use icicle_core::traits::Arithmetic;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;
    use icicle_core::{msm, msm::MSMConfig};

    use crate::{
        kzg::{DegreeCheckVK, PK, VK},
        utils::srs::unsafe_setup_from_tau,
        utils::evaluate_all_lagrange_coefficients
    };

    use super::{
        structs::{Instance, Witness},
        InnerProduct,
    };

    use icicle_core::ntt::{NTTInitDomainConfig, initialize_domain, release_domain};

    const N: usize = 32;
    const LOG_N: usize = 5;

    #[test]
    fn test_ipa() {
        let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

        let domain = get_root_of_unity::<Bn254ScalarField>(N as u64);

        let tau = Bn254ScalarField::from_u32(100u32);
        let lb_at_tau = evaluate_all_lagrange_coefficients::<Bn254CurveCfg>(domain, tau, N);

        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;
        
        let mut x_g2_affine = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&x_g2, &mut x_g2_affine);
        
        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2_affine);

        // we will be checking for R <= n - 2
        let shift_factor = srs.len() - 1 - (N - 2);

        let tau_pow_shift_projective = Bn254G2CurveCfg::mul_scalar(Bn254G2CurveCfg::get_generator(), tau.pow(shift_factor as usize));
        let mut tau_pow_shift_affine = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&tau_pow_shift_projective, &mut tau_pow_shift_affine);
        
        let mut degree_check_vk_map: BTreeMap<usize, Affine<Bn254G2CurveCfg>> = BTreeMap::new();
        degree_check_vk_map.insert(shift_factor, tau_pow_shift_affine);
        let degree_check_vk = DegreeCheckVK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> {
            pk_max_degree: srs.len() - 1,
            shifts: degree_check_vk_map,
            e: PhantomData,
        };
        
        let generator = Bn254CurveCfg::get_generator();
        let a: Vec<Bn254ScalarField> = ScalarCfg::generate_random(N);

        let mut lagrange_basis: Vec<Affine<Bn254CurveCfg>> = Vec::with_capacity(lb_at_tau.len());
        for li in lb_at_tau {
            let point = Bn254CurveCfg::mul_scalar(generator, li);
            lagrange_basis.push(point.into());
        }

        let b = Bn254CurveCfg::generate_random_affine_points(N);
        let h_base_projective = Bn254CurveCfg::mul_scalar(generator, ScalarCfg::generate_random(1)[0]);
        let mut h_base_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&h_base_projective, &mut h_base_affine);

        let cfg = MSMConfig::default();
        let mut ac_projective = vec![Projective::<Bn254CurveCfg>::zero(); 1];
        
        msm::msm(
            HostSlice::from_slice(&a),
            HostSlice::from_slice(&lagrange_basis),
            &cfg,
            HostSlice::from_mut_slice(&mut ac_projective),
        )
        .unwrap();
        
        let mut ac_affine = Affine::<Bn254CurveCfg>::zero();
        Bn254CurveCfg::to_affine(&ac_projective[0], &mut ac_affine);
        
        let mut c: Vec<Affine::<Bn254CurveCfg>> = Vec::with_capacity(b.len());

        for i in 0..b.len() {
            let bi = b[i];
            let ai = a[i];

            let c_projective = Bn254CurveCfg::mul_scalar(bi.to_projective(), ai);
            let mut c_affine = Affine::<Bn254CurveCfg>::zero();
            Bn254CurveCfg::to_affine(&c_projective, &mut c_affine);
            c.push(c_affine);
        }

        let instance = Instance::<N, Bn254CurveCfg> {
            ac: ac_affine,
            b: b.try_into().unwrap(),
            h_base: h_base_affine,
            c: c.try_into().unwrap(),
        };

        let witness = Witness::<N, Bn254CurveCfg> {
            a: a.try_into().unwrap(),
        };

        let proof = InnerProduct::<N, LOG_N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::prove(&instance, &witness, &pk, 0usize).unwrap();
        let res = InnerProduct::<N, LOG_N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&instance, &proof, &vk, &degree_check_vk, 0usize);
        
        release_domain::<Bn254ScalarField>().unwrap();

        assert!(res.is_ok());
    }
}