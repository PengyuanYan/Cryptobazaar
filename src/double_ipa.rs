use icicle_core::curve::{Curve,Affine,Projective};
use icicle_core::traits::FieldImpl;
use icicle_core::pairing::Pairing;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_runtime::memory::HostSlice;
use icicle_core::ntt::{NTTDomain, NTT};
use icicle_core::traits::GenerateRandom;
use icicle_core::traits::Arithmetic;

use icicle_core::{msm, msm::MSMConfig};

use crate::double_pedersen_schnorr::{
    structs::{Instance as PSInstance, Witness as PSWitness},
    Argument as PSArgument,
};
use crate::fold_lagrange::{structs::Instance as LFInstance, Argument as LFArgument};
use crate::kzg::{PK as KzgPk, VK as KzgVk};
use crate::utils::folding::{compute_folding_coeffs, AffFold, FFold, Fold};

use self::{
    structs::{Instance, Proof, Witness},
    tr::Transcript,
};

pub mod structs;
mod tr;

#[derive(Debug)]
pub enum Error {
    FoldShapeMismatch { round: usize, a_left: usize, b1_right: usize, b2_right: usize,
                                      a_right: usize, b1_left: usize, b2_left: usize
                      },
    FoldWrongLength { a_len: usize, b1_len: usize, b2_len: usize },
    SanityMismatch { position: usize },
}

pub struct DoubleInnerProduct<const N: usize, const LOG_N: usize, C1, C2, F, U>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    _e: PhantomData<(C1, C2, F, U)>,
}

impl<const N: usize, const LOG_N: usize, C1, C2, F, U> DoubleInnerProduct<N, LOG_N, C1, C2, F, U>
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
    ) -> Result<Proof<LOG_N, C1>, Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: GenerateRandom<<C1 as Curve>::ScalarField>,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        let mut acc_blinders_1 = C1::ScalarField::zero();
        let mut acc_blinders_2 = C1::ScalarField::zero();

        let mut l_1_msgs = Vec::<Affine::<C1>>::with_capacity(LOG_N);
        let mut r_1_msgs = Vec::<Affine::<C1>>::with_capacity(LOG_N);
        let mut l_2_msgs = Vec::<Affine::<C1>>::with_capacity(LOG_N);
        let mut r_2_msgs = Vec::<Affine::<C1>>::with_capacity(LOG_N);

        let mut alpha_invs = Vec::with_capacity(LOG_N);

        let mut tr = Transcript::<N, LOG_N, C1>::new(b"double-ipa");
        tr.send_instance(instance);

        let r = tr.get_r();
        let r_pows: Vec<_> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * r)).take(N).collect();

        let mut c1_fold_projective = instance.ac.clone().to_projective();

        let cfg = MSMConfig::default();
        let mut c2_fold_projective_v = vec![Projective::<C1>::zero(); 1];
        msm::msm(
            HostSlice::from_slice(&r_pows),
            HostSlice::from_slice(&instance.c),
            &cfg,
            HostSlice::from_mut_slice(&mut c2_fold_projective_v),
        )
        .unwrap();

        let mut c2_fold_projective = c2_fold_projective_v[0];

        let mut a_folded = witness.a.clone().to_vec();

        let mut b1_folded = instance.lagrange_basis.clone().to_vec();

        let mut b2_folded: Vec<Affine::<C1>> = Vec::with_capacity(instance.b.len());
        for i in 0..instance.b.len() {
            let bi = instance.b[i];
            let ri = r_pows[i];

            let result_projective = C1::mul_scalar(bi.to_projective(), ri);

            let mut result_affine = Affine::<C1>::zero();
            C1::to_affine(&result_projective, &mut result_affine);

            b2_folded.push(result_affine);
        }
        
        for round in 0..LOG_N as usize {
            let a_left = &a_folded[..a_folded.len() / 2];
            let a_right = &a_folded[a_folded.len() / 2..];

            let b1_left = &b1_folded[..b1_folded.len() / 2];
            let b1_right = &b1_folded[b1_folded.len() / 2..];

            let b2_left = &b2_folded[..b2_folded.len() / 2];
            let b2_right = &b2_folded[b2_folded.len() / 2..];

            assert_eq!(a_left.len(), b1_right.len());
            assert_eq!(a_left.len(), b2_right.len());
            assert_eq!(a_right.len(), b1_left.len());
            assert_eq!(a_right.len(), b2_left.len());

            if a_left.len() != b1_right.len() ||
               a_left.len() != b2_right.len() ||

               a_right.len() != b1_left.len() ||
               a_right.len() != b2_left.len()
            {
                return Err(Error::FoldShapeMismatch {
                    round: round,
                    a_left: a_left.len(),
                    b1_right: b1_right.len(),
                    b2_right: b2_right.len(),

                    a_right: a_right.len(),
                    b1_left: b1_left.len(),
                    b2_left: b2_left.len(),
                });
            }

            let l_1 = if a_left.len() == 1 { 
                C1::mul_scalar(b1_right[0].to_projective(), a_left[0])
            } else {
                let mut buf = [Projective::<C1>::zero()];
                msm::msm(
                    HostSlice::from_slice(a_left),
                    HostSlice::from_slice(b1_right),
                    &cfg,
                    HostSlice::from_mut_slice(&mut buf),
                )
                .unwrap();
                buf[0]
            };

            let r_1 = if a_right.len() == 1 {
                C1::mul_scalar(b1_left[0].to_projective(), a_right[0])
            } else {
                let mut buf = [Projective::<C1>::zero()];
                msm::msm(
                    HostSlice::from_slice(a_right),
                    HostSlice::from_slice(b1_left),
                    &cfg,
                    HostSlice::from_mut_slice(&mut buf),
                )
                .unwrap();
                buf[0]
            };

            let l_2 = if a_left.len() == 1 { 
                C1::mul_scalar(b2_right[0].to_projective(), a_left[0])
            } else {
                let mut buf = [Projective::<C1>::zero()];
                msm::msm(
                    HostSlice::from_slice(a_left),
                    HostSlice::from_slice(b2_right),
                    &cfg,
                    HostSlice::from_mut_slice(&mut buf),
                )
                .unwrap();
                buf[0]
            };

            let r_2 = if a_right.len() == 1 {
                C1::mul_scalar(b2_left[0].to_projective(), a_right[0])
            } else {
                let mut buf = [Projective::<C1>::zero()];
                msm::msm(
                    HostSlice::from_slice(a_right),
                    HostSlice::from_slice(b2_left),
                    &cfg,
                    HostSlice::from_mut_slice(&mut buf),
                )
                .unwrap();
                buf[0]
            };

            let blinders = <<C1::ScalarField as FieldImpl>::Config as GenerateRandom<C1::ScalarField>>::generate_random(4);
            let blinder_l1 = blinders[0];
            let blinder_r1 = blinders[1];

            let blinder_l2 = blinders[2];
            let blinder_r2 = blinders[3];

            let l_1_projective = l_1 + C1::mul_scalar(instance.h_base.to_projective(), blinder_l1);
            let mut l_1_affine = Affine::<C1>::zero();
            C1::to_affine(&l_1_projective, &mut l_1_affine);
            l_1_msgs.push(l_1_affine);

            let r_1_projective = r_1 + C1::mul_scalar(instance.h_base.to_projective(), blinder_r1);
            let mut r_1_affine = Affine::<C1>::zero();
            C1::to_affine(&r_1_projective, &mut r_1_affine);
            r_1_msgs.push(r_1_affine);

            let l_2_projective = l_2 + C1::mul_scalar(instance.h_base.to_projective(), blinder_l2);
            let mut l_2_affine = Affine::<C1>::zero();
            C1::to_affine(&l_2_projective, &mut l_2_affine);
            l_2_msgs.push(l_2_affine);

            let r_2_projective = r_2 + C1::mul_scalar(instance.h_base.to_projective(), blinder_r2);
            let mut r_2_affine = Affine::<C1>::zero();
            C1::to_affine(&r_2_projective, &mut r_2_affine);
            r_2_msgs.push(r_2_affine);

            tr.send_ls_rs(&l_1_affine, &r_1_affine, &l_2_affine, &r_2_affine);
            let alpha = tr.get_alpha_i();
            let alpha_inv = alpha.inv();
            alpha_invs.push(alpha_inv);

            acc_blinders_1 = acc_blinders_1 + alpha_inv * blinder_l1 + alpha * blinder_r1;
            acc_blinders_2 = acc_blinders_2 + alpha_inv * blinder_l2 + alpha * blinder_r2;

            // fold vectors
            a_folded = FFold::<C1>::fold_vec(&a_folded, alpha).unwrap();
            b1_folded = AffFold::<C1>::fold_vec(&b1_folded, alpha_inv).unwrap();
            b2_folded = AffFold::<C1>::fold_vec(&b2_folded, alpha_inv).unwrap();

            // derive new cm
            c1_fold_projective = C1::mul_scalar(l_1_projective, alpha_inv) + c1_fold_projective + C1::mul_scalar(r_1_projective, alpha);
            c2_fold_projective = C1::mul_scalar(l_2_projective, alpha_inv) + c2_fold_projective + C1::mul_scalar(r_2_projective, alpha);
        }
        
        // sanity
        if a_folded.len() != 1 || 
           b1_folded.len() != 1 ||
           b2_folded.len() != 1
        {
            return Err(Error::FoldWrongLength {
                a_len: a_folded.len(),
                b1_len: b1_folded.len(),
                b2_len: b2_folded.len(),
            });
        }

        let lhs_1 = C1::mul_scalar(b1_folded[0].to_projective(), a_folded[0])
                  + C1::mul_scalar(instance.h_base.to_projective(), acc_blinders_1);
        let lhs_2 = C1::mul_scalar(b2_folded[0].to_projective(), a_folded[0])
                  + C1::mul_scalar(instance.h_base.to_projective(), acc_blinders_2);

        if lhs_1 != c1_fold_projective {
            return Err(Error::SanityMismatch {
                position: 1usize
            });
        }

        if lhs_2 != c2_fold_projective {
            return Err(Error::SanityMismatch {
                position: 2usize
            });
        }
        
        let mut c1_fold_affine = Affine::<C1>::zero();
        C1::to_affine(&c1_fold_projective, &mut c1_fold_affine);
        let mut c2_fold_affine = Affine::<C1>::zero();
        C1::to_affine(&c2_fold_projective, &mut c2_fold_affine);

        let ps_instance = PSInstance::<C1> {
            q_base: b1_folded[0],
            p_base: b2_folded[0],
            h_base: instance.h_base,
            x_1: c1_fold_affine,
            x_2: c2_fold_affine,
        };

        let ps_witness = PSWitness {
            a: a_folded[0],
            r_1: acc_blinders_1,
            r_2: acc_blinders_2,
        };
        
        let ps_proof = PSArgument::prove(&ps_instance, &ps_witness);
        
        let lf_instance = LFInstance::<N, LOG_N, C1> {
            lb_commitments: instance.lagrange_basis,
            challenges: alpha_invs.try_into().unwrap(),
        };
        
        let lf_proof = LFArgument::<N, LOG_N, C1, C2, F, U>::prove(pk, &lf_instance);
        
        Ok(Proof::<LOG_N, C1> {
            l_1: l_1_msgs.try_into().unwrap(),
            r_1: r_1_msgs.try_into().unwrap(),
            l_2: l_2_msgs.try_into().unwrap(),
            r_2: r_2_msgs.try_into().unwrap(),
            lf_proof,
            ps_proof,
        })
    }

    pub fn verify(
        instance: &Instance<N, C1>,
        proof: &Proof<LOG_N, C1>,
        vk: &KzgVk<C1, C2, F>
    )
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>,
    {
        let mut tr = Transcript::<N, LOG_N, _>::new(b"double-ipa");
        tr.send_instance(instance);

        let r = tr.get_r();
        let r_pows: Vec<_> = std::iter::successors(Some(C1::ScalarField::one()), |p| Some(*p * r)).take(N).collect();

        let mut c1_fold_projective = instance.ac.clone().to_projective();

        let cfg = MSMConfig::default();
        let mut c2_fold_projective_v = vec![Projective::<C1>::zero(); 1];
        msm::msm(
            HostSlice::from_slice(&r_pows),
            HostSlice::from_slice(&instance.c),
            &cfg,
            HostSlice::from_mut_slice(&mut c2_fold_projective_v),
        )
        .unwrap();

        let mut c2_fold_projective = c2_fold_projective_v[0];
        
        let mut b2_rescaled: Vec<Affine::<C1>> = Vec::with_capacity(instance.b.len());

        for i in 0..instance.b.len() {
            let bi = instance.b[i];
            let ri = r_pows[i];

            let result_projective = C1::mul_scalar(bi.to_projective(), ri);

            let mut result_affine = Affine::<C1>::zero();
            C1::to_affine(&result_projective, &mut result_affine);

            b2_rescaled.push(result_affine);
        }

        let mut alpha_invs = Vec::with_capacity(LOG_N);

        for i in 0..LOG_N as usize {
            tr.send_ls_rs(&proof.l_1[i], &proof.r_1[i], &proof.l_2[i], &proof.r_2[i]);

            let alpha = tr.get_alpha_i();
            let alpha_inv = alpha.inv();
            alpha_invs.push(alpha_inv);

            c1_fold_projective = C1::mul_scalar(proof.l_1[i].to_projective(), alpha_inv) + c1_fold_projective + C1::mul_scalar(proof.r_1[i].to_projective(), alpha);
            c2_fold_projective = C1::mul_scalar(proof.l_2[i].to_projective(), alpha_inv) + c2_fold_projective + C1::mul_scalar(proof.r_2[i].to_projective(), alpha);
        }

        // now fold the basis b and check pedersen openings
        let fold_coeffs = compute_folding_coeffs::<C1>(&alpha_invs);
        let mut b2_folded_projective_v = vec![Projective::<C1>::zero(); 1];
        msm::msm(
            HostSlice::from_slice(&fold_coeffs),
            HostSlice::from_slice(&b2_rescaled),
            &cfg,
            HostSlice::from_mut_slice(&mut b2_folded_projective_v),
        )
        .unwrap();

        let b2_folded_projective = b2_folded_projective_v[0];
        let mut b2_folded_affine = Affine::<C1>::zero();
        C1::to_affine(&b2_folded_projective, &mut b2_folded_affine);

        let lf_instance = LFInstance::<N, LOG_N, C1> {
            lb_commitments: instance.lagrange_basis,
            challenges: alpha_invs.try_into().unwrap(),
        };
        
        let lf_folding_res = LFArgument::<N, LOG_N, C1, C2, F, U>::verify(&proof.lf_proof, &lf_instance, vk);
        assert!(lf_folding_res.is_ok());

        let mut c1_fold_affine = Affine::<C1>::zero();
        C1::to_affine(&c1_fold_projective, &mut c1_fold_affine);
        let mut c2_fold_affine = Affine::<C1>::zero();
        C1::to_affine(&c2_fold_projective, &mut c2_fold_affine);

        let ps_instance = PSInstance::<C1> {
            q_base: proof.lf_proof.p_cm,
            p_base: b2_folded_affine,
            h_base: instance.h_base,
            x_1: c1_fold_affine,
            x_2: c2_fold_affine,
        };

        let res = PSArgument::verify(&ps_instance, &proof.ps_proof);
        assert!(res.is_ok());
    }
}

#[cfg(test)]
mod double_ipa_tests {
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::curve::{Curve,Affine,Projective};
    use icicle_core::traits::FieldImpl;
    use std::marker::PhantomData;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;
    use icicle_core::ntt::get_root_of_unity;
    use icicle_runtime::memory::HostSlice;
    use icicle_bn254::curve::ScalarCfg;
    use icicle_core::traits::GenerateRandom;
    use icicle_core::{msm, msm::MSMConfig};

    use crate::{
        kzg::{PK, VK},
        utils::srs::unsafe_setup_from_tau,
        utils::evaluate_all_lagrange_coefficients,
    };

    use super::{
        structs::{Instance, Witness},
        DoubleInnerProduct,
    };
    
    use icicle_core::ntt::{NTTInitDomainConfig, initialize_domain, release_domain};

    const N: usize = 32;
    const LOG_N: usize = 5;

    #[test]
    fn test_double_ipa() {
        let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

        let tau = Bn254ScalarField::from_u32(100u32);
        let domain = get_root_of_unity::<Bn254ScalarField>(N as u64);
        let lb_at_tau = evaluate_all_lagrange_coefficients::<Bn254CurveCfg>(domain, tau, N);

        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let mut x_g2_affine = Affine::<Bn254G2CurveCfg>::zero();
        Bn254G2CurveCfg::to_affine(&x_g2, &mut x_g2_affine);

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2_affine);

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
            lagrange_basis: lagrange_basis.try_into().unwrap(),
            b: b.try_into().unwrap(),
            h_base: h_base_affine,
            c: c.try_into().unwrap(),
        };

        let witness = Witness::<N, Bn254CurveCfg> {
            a: a.try_into().unwrap(),
        };
        
        let proof =
            DoubleInnerProduct::<N, LOG_N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::prove(&instance, &witness, &pk).unwrap();
        
        DoubleInnerProduct::<N, LOG_N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&instance, &proof, &vk);

        release_domain::<Bn254ScalarField>().unwrap();
    }
}