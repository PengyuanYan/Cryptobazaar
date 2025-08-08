use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{ntt, NTTInitDomainConfig, NTTConfig, NTTDir, get_root_of_unity, NTTDomain, NTT};
use icicle_core::traits::Arithmetic;
use icicle_runtime::{Device};//,runtime};

use icicle_core::curve::{Curve, Projective, Affine};
use icicle_core::{msm, msm::MSMConfig};
use icicle_runtime::memory::{DeviceVec, HostSlice};

pub mod folding;
pub mod srs;

pub fn is_pow_2(x: usize) -> bool {
    (x & (x - 1)) == 0
}

pub fn evaluate_vanishing_over_extended_coset<C: Curve>(n: usize, k: usize) -> Vec<C::ScalarField> 
where 
    <C as Curve>::ScalarField: Arithmetic,
    <<C as Curve>::ScalarField as FieldImpl>::Config: NTTDomain<<C as Curve>::ScalarField>,
{
    assert!(is_pow_2(n));
    assert!(is_pow_2(k));

    let p = C::ScalarField::from_u32(5u32);
    let domain_kn = get_root_of_unity::<C::ScalarField>((k * n).try_into().unwrap());

    let coset_generator_pow_n = p.pow(n.try_into().unwrap());
    let wi = domain_kn;

    let mut modulus_zh_coset_evals = Vec::with_capacity(k);

    for i in 0usize..k {
        let zhi = coset_generator_pow_n * wi.pow((i * n).try_into().unwrap()) - C::ScalarField::one();
        modulus_zh_coset_evals.push(zhi);
    }

    modulus_zh_coset_evals
}

pub fn get_coeffs_of_poly<P>(poly: &P) -> Vec<P::Field>
where
    P: UnivariatePolynomial,
    P::Field: FieldImpl,
{
    let n: usize = poly.get_nof_coeffs().try_into().expect("Host's archtecture is smaller than 64-bit");

    let mut coeffs = vec![P::Field::zero(); n];

    poly.copy_coeffs(0, HostSlice::from_mut_slice(&mut coeffs));

    coeffs
}

// this function is reimplemented based on the source code of ark library
pub fn evaluate_all_lagrange_coefficients<C: Curve>(domain: C::ScalarField, tau: C::ScalarField, n: usize) -> Vec<C::ScalarField>
where
    <C as Curve>::ScalarField: Arithmetic,
{
    let size = n;
    let z_h_at_tau = tau.pow(n) - C::ScalarField::one();
    let offset = C::ScalarField::one();
    let group_gen = domain;

    assert_eq!(domain.pow(n), C::ScalarField::one());

    if z_h_at_tau == C::ScalarField::zero() {
        let mut u = vec![C::ScalarField::zero(); size];
        let mut omega_i = offset;
        for u_i in u.iter_mut().take(size) {
            if omega_i == tau {
                *u_i = C::ScalarField::one();
                break;
            }
            omega_i = omega_i * group_gen;
        }
        
        u

    } else {
        let group_gen_inv = domain.inv();

        // v_0_inv = m * h^(m-1)
        let v_0_inv = C::ScalarField::from_u32(n as u32) * offset.pow(size - 1);
        let mut l_i = z_h_at_tau.inv() * v_0_inv;
        let mut negative_cur_elem = C::ScalarField::zero() - offset;
        let mut lagrange_coefficients_inverse = vec![C::ScalarField::zero(); size];
        for coeff in &mut lagrange_coefficients_inverse {
            let r_i = tau + negative_cur_elem;
            *coeff = l_i * r_i;
            // Increment l_i and negative_cur_elem
            l_i = l_i * group_gen_inv;
            negative_cur_elem = negative_cur_elem * group_gen;
        }

        // Invert the lagrange coefficients inverse, to get the actual coefficients,
        // and return these
        for i in 0..lagrange_coefficients_inverse.len() {
            lagrange_coefficients_inverse[i] = lagrange_coefficients_inverse[i].inv();
        }

        lagrange_coefficients_inverse
    }
}

pub fn load_backend() {
    //runtime::load_backend_from_env_or_default().unwrap();
    let device_gpu = Device::new("CUDA", 0);
    let is_cuda_device_available = icicle_runtime::is_device_available(&device_gpu);
    if is_cuda_device_available {
        println!("gpu");
        icicle_runtime::set_device(&device_gpu).unwrap();
        let device = icicle_runtime::get_active_device().unwrap();
        if device.get_device_type() == "CUDA" {
            println!("gpu correct");
        } else {
            println!("gpu false");
        }
    } else {
        println!("cpu");
        let device_cpu = Device::new("CPU", 0);
        icicle_runtime::set_device(&device_cpu).unwrap();
        let device = icicle_runtime::get_active_device().unwrap();
        if device.get_device_type() == "CPU" {
            println!("cpu correct");
        } else {
            println!("cpu false");
        }
    }
}

pub fn get_device_is_cpu_or_gpu() -> usize {
    let mut cpu_or_gpu: usize = 0;
    let device = icicle_runtime::get_active_device().unwrap();
    if device.get_device_type() == "CUDA" {
        println!("gpu correct");
        cpu_or_gpu = 1;
    } else if device.get_device_type() == "CPU"{
        println!("cpu correct");
    } else {
        panic!("cant find gpu or cpu");
    }

    return cpu_or_gpu;
}

pub fn msm_gpu<C: Curve + icicle_core::msm::MSM<C>>(scalars: &[C::ScalarField], bases: &[Affine::<C>], cfg: &MSMConfig) -> Projective::<C> {
    assert_eq!(scalars.len(), bases.len());
    
    let cfg = MSMConfig::default();
    let n: usize = scalars.len();
    let mut output_cpu = vec![Projective::<C>::zero(); 1];

    let mut input_scalars_gpu =
        DeviceVec::<C::ScalarField>::device_malloc(n).expect("Failed to allocate device memory for scalars input");
    let mut input_bases_gpu =
        DeviceVec::<Affine::<C>>::device_malloc(n).expect("Failed to allocate device memory for bases input");
    let mut output_gpu =
        DeviceVec::<Projective::<C>>::device_malloc(1).expect("Failed to allocate device memory for msm output");
    
    input_scalars_gpu
        .copy_from_host(HostSlice::from_slice(&scalars))
        .expect("Failed to copy scalars data to GPU");
    input_bases_gpu
        .copy_from_host(HostSlice::from_slice(&bases))
        .expect("Failed to copy bases data to GPU");

    msm::msm(
        &input_scalars_gpu[..],
        &input_bases_gpu[..],
        &cfg,
        &mut output_gpu[..],
    )
    .unwrap();

    output_gpu
        .copy_to_host(HostSlice::from_mut_slice(&mut output_cpu))
        .expect("Failed to copy result data back to CPU");
    
    output_cpu[0]
}

pub fn ntt_gpu<C: Curve>(scalars: &[C::ScalarField], cfg: &NTTConfig<C::ScalarField>, direction: NTTDir) -> Vec<C::ScalarField>
where
    <<C as Curve>::ScalarField as FieldImpl>::Config: NTT<<C as Curve>::ScalarField, <C as Curve>::ScalarField>
{
    let n = scalars.len();

    let mut input_gpu =
        DeviceVec::<C::ScalarField>::device_malloc(n).expect("Failed to allocate device memory for ntt input");
    let mut output_gpu =
        DeviceVec::<C::ScalarField>::device_malloc(n).expect("Failed to allocate device memory for ntt output");
    input_gpu
        .copy_from_host(HostSlice::from_slice(&scalars))
        .expect("Failed to copy ntt input data to GPU");
    
    let mut coeffs = vec![C::ScalarField::zero(); n];
    ntt(
        &input_gpu[..],
        direction,
        cfg,
        &mut output_gpu[..],
    )
    .unwrap();

    output_gpu
        .copy_to_host(HostSlice::from_mut_slice(&mut coeffs))
        .expect("Failed to copy ntt result data back to CPU");

    coeffs
}

pub fn my_ntt<C: Curve>(scalars: &[C::ScalarField], coset_gen: C::ScalarField, direction: NTTDir, cpu_or_gpu: usize) -> Vec<C::ScalarField>
where
    <<C as Curve>::ScalarField as FieldImpl>::Config: NTT<<C as Curve>::ScalarField, <C as Curve>::ScalarField>
{
    let mut cfg = NTTConfig::<C::ScalarField>::default();

    if coset_gen != C::ScalarField::one() {
        cfg.coset_gen = coset_gen;
        println!("not one");
    } else {
        println!("is one");
    }

    let result_coeffs = if cpu_or_gpu == 0 {
        let mut output_scalars = vec![C::ScalarField::zero(); scalars.len()];
            ntt(
                HostSlice::from_slice(scalars),
                direction,
                &cfg,
                HostSlice::from_mut_slice(&mut output_scalars),
            )
            .unwrap();

            output_scalars
        } else {
            ntt_gpu::<C>(scalars, &cfg, direction)
    };

    result_coeffs
}

pub fn my_msm<C: Curve + icicle_core::msm::MSM<C>>(scalars: &[C::ScalarField], bases: &[Affine::<C>], cpu_or_gpu: usize) -> Projective::<C> {
    let cfg = MSMConfig::default();

    let msm_result = if cpu_or_gpu == 0 {
            let mut msm_result_v = [Projective::<C>::zero()];
                msm::msm(
                    HostSlice::from_slice(scalars),
                    HostSlice::from_slice(bases),
                    &cfg,
                    HostSlice::from_mut_slice(&mut msm_result_v),
                )
                .unwrap();

                msm_result_v[0]
            } else {
                msm_gpu(scalars, bases, &cfg)
    };

    msm_result
}