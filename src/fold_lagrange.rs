/*
    1. P sends [p] claimed to be correct folding of lagrange basis commitments given array of challenges {ch}
    2. V sends ß
    3. P
        1. computes a(X) = ∑ß^iLi(X)
        2. sends [a] and π (well formation proof of a)

    4: P, V
        Run univariate sumcheck to prove that ∑p_i * a_i = fold(ß)
*/

use icicle_core::curve::{Curve, Affine, Projective};
use icicle_core::traits::FieldImpl;
use icicle_core::pairing::Pairing;
use std::marker::PhantomData;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_runtime::memory::HostSlice;
use icicle_core::ntt::{ntt, NTTDomain, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, initialize_domain, NTT, release_domain};
use icicle_core::traits::Arithmetic;
use icicle_core::{msm, msm::MSMConfig};

use self::structs::{Error, Proof};
use self::{structs::Instance, tr::Transcript};
use crate::acc::{
    structs::{Instance as AccInstance, Witness as AccWitness},
    Argument as AccArgument,
};

use crate::univariate_sumcheck::structs::{Instance as UVInstance, Witness as UVWitness};
use crate::univariate_sumcheck::UnivariateSumcheck;
use crate::{
    kzg::{Kzg, PK as KzgPk, VK as KzgVk},
    utils::folding::compute_folding_coeffs,
};

pub mod structs;
mod tr;

pub struct Argument<const N: usize, const LOG_N: usize, C1, C2, F, U>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    _e: PhantomData<(C1, C2, F, U)>,
}

impl<const N: usize, const LOG_N: usize, C1, C2, F, U> Argument<N, LOG_N, C1, C2, F, U>
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = <C1 as Curve>::ScalarField>,
{
    pub fn prove(pk: &KzgPk<C1, C2, F>, instance: &Instance<N, LOG_N, C1>) -> Proof<C1>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
    {
        let mut tr = Transcript::<C1>::new(b"fold-lagrange");

        let folding_coeffs = compute_folding_coeffs::<C1>(&instance.challenges);
        
        let cfg_msm = MSMConfig::default();
        let mut p_cm_buf = [Projective::<C1>::zero()];
        msm::msm(
            HostSlice::from_slice(&folding_coeffs),
            HostSlice::from_slice(&instance.lb_commitments),
            &cfg_msm,
            HostSlice::from_mut_slice(&mut p_cm_buf),
        )
        .unwrap();
        let p_cm_projective = p_cm_buf[0];

        let mut p_cm_affine = Affine::<C1>::zero();
        C1::to_affine(&p_cm_projective, &mut p_cm_affine);

        tr.send_p(&p_cm_affine);
        let beta = tr.get_beta();

        let mut beta_powers: Vec<C1::ScalarField> = Vec::with_capacity(N);
        let mut power = C1::ScalarField::one();
        for _ in 0..N {
            beta_powers.push(power);
            power = power * beta;
        }

        let mut sum = C1::ScalarField::zero();
        for i in 0..folding_coeffs.len() {
            sum = folding_coeffs[i] * beta_powers[i] + sum;
        }
        
        //let domain = get_root_of_unity::<C1::ScalarField>((N).try_into().unwrap());
        let cfg_ntt = NTTConfig::<C1::ScalarField>::default();
        //initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();
        let mut p_coeffs = vec![C1::ScalarField::zero(); folding_coeffs.len()];
        ntt(
            HostSlice::from_slice(&folding_coeffs),
            NTTDir::kInverse,
            &cfg_ntt,
            HostSlice::from_mut_slice(&mut p_coeffs),
        )
        .unwrap();
        
        let p = U::from_coeffs(HostSlice::from_slice(&p_coeffs), p_coeffs.len());
        
        let mut acc_coeffs = vec![C1::ScalarField::zero(); beta_powers.len()];
        ntt(
            HostSlice::from_slice(&beta_powers),
            NTTDir::kInverse,
            &cfg_ntt,
            HostSlice::from_mut_slice(&mut acc_coeffs),
        )
        .unwrap();
        
        //release_domain::<C1::ScalarField>().unwrap();

        let acc = U::from_coeffs(HostSlice::from_slice(&acc_coeffs), acc_coeffs.len());

        let acc_cm = Kzg::commit(pk, &acc).unwrap();

        /*** Sumcheck */
        let uv_instance = UVInstance::<C1> {
            n: N,
            a_cm: p_cm_affine.clone(),
            b_cm: acc_cm.clone(),
            sum,
        };

        let uv_witness = UVWitness::<C1, U> {
            a_poly: p,
            b_poly: acc.clone(),
            e: PhantomData,
        };

        let sumcheck_proof = UnivariateSumcheck::<C1, C2, F, U>::prove(&uv_witness, &uv_instance, pk);
        /*** Sumcheck */

        /*** Acc */
        let acc_instance = AccInstance::<C1> {
            n: N,
            mu: beta,
            acc_cm,
        };

        let acc_witness = AccWitness::<C1, U> { acc, e: PhantomData, };

        let acc_proof = AccArgument::<C1, C2, F>::prove(&acc_instance, &acc_witness, pk);
        /*** Acc */

        Proof {
            p_cm: p_cm_affine,
            acc_cm,
            acc_proof,
            sumcheck_proof,
        }
    }

    pub fn verify(
        proof: &Proof<C1>,
        instance: &Instance<N, LOG_N, C1>,
        vk: &KzgVk<C1, C2, F>,
    ) -> Result<(), Error>
    where
        <C1 as Curve>::ScalarField: Arithmetic,
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C1 as Curve>::ScalarField>
    {
        let mut tr = Transcript::<C1>::new(b"fold-lagrange");
        let folding_coeffs = compute_folding_coeffs::<C1>(&instance.challenges);

        tr.send_p(&proof.p_cm);
        let beta = tr.get_beta();

        let acc_instance = AccInstance::<C1> {
            n: N,
            mu: beta,
            acc_cm: proof.acc_cm,
        };
        
        let mut beta_powers: Vec<C1::ScalarField> = Vec::with_capacity(N);
        let mut power = C1::ScalarField::one();
        for _ in 0..N {
            beta_powers.push(power);
            power = power * beta;
        }
        
        let mut sum = C1::ScalarField::zero();
        for i in 0..folding_coeffs.len() {
            sum = folding_coeffs[i] * beta_powers[i] + sum;
        }
        
        let uv_instance = UVInstance::<C1> {
            n: N,
            a_cm: proof.p_cm,
            b_cm: proof.acc_cm,
            sum,
        };

        UnivariateSumcheck::<C1, C2, F, U>::verify(&proof.sumcheck_proof, &uv_instance, vk)?;
        AccArgument::<C1, C2, F>::verify(&acc_instance, &proof.acc_proof, vk)?;
        Ok(())
    }
}

#[cfg(test)]
mod lagrange_fold_tests {
    use icicle_bn254::curve::{CurveCfg as Bn254CurveCfg, G2CurveCfg as Bn254G2CurveCfg, G1Affine};
    use icicle_bn254::pairing::PairingTargetField as Bn254PairingFieldImpl;
    use icicle_core::curve::Curve;
    use icicle_bn254::polynomials::DensePolynomial as Bn254Poly;

    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;

    use icicle_core::ntt::get_root_of_unity;
    use std::marker::PhantomData;

    use crate::{
        kzg::{PK, VK},
        utils::srs::unsafe_setup_from_tau,
        utils::evaluate_all_lagrange_coefficients,
    };

    use super::{structs::Instance, Argument};

    use icicle_core::ntt::{NTTInitDomainConfig, initialize_domain, release_domain};

    const N: usize = 16;
    const LOG_N: usize = 4;

    #[test]
    fn test_folding() {
        let domain = get_root_of_unity::<Bn254ScalarField>((N * N).try_into().unwrap());
        initialize_domain(domain, &NTTInitDomainConfig::default()).unwrap();

        let g = Bn254CurveCfg::get_generator();
        let tau = Bn254ScalarField::from_u32(17u32);
        let srs = unsafe_setup_from_tau::<Bn254CurveCfg>(N - 1, tau);
        let x_g2 = Bn254G2CurveCfg::get_generator() * tau;

        let pk = PK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl> { srs: srs.clone(), e: PhantomData, };
        let vk = VK::<Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl>::new(x_g2.into());
        
        let domain = get_root_of_unity::<Bn254ScalarField>((N).try_into().unwrap());
        let lb_at_tau = evaluate_all_lagrange_coefficients::<Bn254CurveCfg>(domain, tau, N);
        
        let mut lb_commitments = Vec::with_capacity(lb_at_tau.len());
        for li in lb_at_tau {
            let projective_base_cm = Bn254CurveCfg::mul_scalar(g, li);
            let mut affine_base_cm: G1Affine = G1Affine::zero();
            Bn254CurveCfg::to_affine(&projective_base_cm, &mut affine_base_cm);
            lb_commitments.push(affine_base_cm);
        }

        let mut chs = Vec::with_capacity(LOG_N);
        for i in 0..LOG_N {
            chs.push(Bn254ScalarField::from_u32((10 + i) as u32));
        }

        let instance = Instance::<N, LOG_N, Bn254CurveCfg> {
            lb_commitments: lb_commitments.try_into().unwrap(),
            challenges: chs.try_into().unwrap(),
        };

        let proof = Argument::<N, LOG_N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::prove(&pk, &instance);
        let res = Argument::<N, LOG_N, Bn254CurveCfg, Bn254G2CurveCfg, Bn254PairingFieldImpl, Bn254Poly>::verify(&proof, &instance, &vk);
        
        release_domain::<Bn254ScalarField>().unwrap();

        assert!(res.is_ok());
    }
}