use icicle_core::curve::{Curve, Affine, Projective};
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;

// Generate a structured reference string (SRS) 
// from a trapdoor value tau.
pub fn unsafe_setup_from_tau<C>(
     max_power: usize,
     tau: C::ScalarField,
) -> Vec<Affine<C>>
where
    C: Curve,
    <C as Curve>::ScalarField: Arithmetic,
{
    let size = max_power + 1;
    
    // Compute [1, τ, τ², …, τ^max_power]
    let mut scalars = Vec::with_capacity(size);
    let mut acc = C::ScalarField::from_u32(1u32);
    for _ in 0..size {
        scalars.push(acc);
        acc = acc * tau;
    }
    
    let projective_g1: Projective<C> = C::get_generator();
    
    // Multiply generator by each τ^i and convert to affine
    let mut srs = Vec::with_capacity(size);
    for i in &scalars {
        let projective_base = C::mul_scalar(projective_g1, *i);
        let mut affine_base = Affine::<C>::zero();
        C::to_affine(&projective_base, &mut affine_base);
        srs.push(affine_base);
    }
    
    srs
}