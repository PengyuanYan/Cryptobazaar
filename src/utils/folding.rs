// Folding utilities used by inner-product arguments
use icicle_core::curve::{Curve, Affine};
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;
use crate::utils::to_affine_batched;

use super::is_pow_2;

#[derive(Debug)]
pub enum FoldError {
    TooShort,
    NotPowTwo,
}

pub trait Fold {
    type Challenge;
    type FoldType;
    fn check(x: &[Self::FoldType]) -> Result<(), FoldError> {
        if x.len() < 2 {
            return Err(FoldError::TooShort);
        }

        if !is_pow_2(x.len()) {
            return Err(FoldError::NotPowTwo);
        }

        Ok(())
    }
    fn fold_vec(
        x: &[Self::FoldType],
        ch: Self::Challenge,
    ) -> Result<Vec<Self::FoldType>, FoldError>;
}

// Scalar folding.
//
// Given vector `[l_0,...,l_{m-1}, r_0,...,r_{m-1}]` and challenge `ch`,
// computes `[l_i + ch * r_i]` for all i.
pub struct FFold<C: Curve>(C);

impl<C: Curve> Fold for FFold<C>
where
    <C as Curve>::ScalarField: Arithmetic,
{
    type Challenge = C::ScalarField;
    type FoldType = C::ScalarField;
    fn fold_vec(
        x: &[Self::FoldType],
        ch: Self::Challenge,
    ) -> Result<Vec<Self::FoldType>, FoldError> {
        Self::check(x)?;
        
        let mid = x.len() / 2;
        let left = &x[..mid];
        let right = &x[mid..];

        let mut folded = Vec::with_capacity(left.len());

        for i in 0..left.len() {
            let l = left[i];
            let r = right[i];
            
            // l_i + ch * r_i
            folded.push(l + r * ch);
        }

        Ok(folded)
    }
}

/// Curve point folding.
///
/// Given points `[L_0,...,L_{m-1}, R_0,...,R_{m-1}]` and challenge `ch`,
/// computes `[L_i + ch * R_i]` for all i.
///
/// The fold is done in projective form for efficiency, then
/// converted back to affine using a batched conversion.
pub struct AffFold<C: Curve>(C);

impl<C: Curve> Fold for AffFold<C> {
    type Challenge = C::ScalarField;
    type FoldType = Affine::<C>;
    fn fold_vec(
        x: &[Self::FoldType],
        ch: Self::Challenge,
    ) -> Result<Vec<Self::FoldType>, FoldError> {
        Self::check(x)?;

        let mid = x.len() / 2;
        let left = &x[..mid];
        let right = &x[mid..];

        let affine_result: Vec<Affine::<C>> = {
            let mut tmp = Vec::with_capacity(left.len());
            for i in 0..left.len() {
                // L_i + ch * R_i
                let result_projective = left[i].to_projective() + C::mul_scalar(right[i].to_projective(), ch);
                tmp.push(result_projective);
            }
            // Batch convert projectives to affines
            to_affine_batched(&tmp)
        };

        Ok(affine_result)
    }
}

// Compute all folding coefficients given a sequence of challenges.
// It is used to compute the \produt (1+\alpha) after the IPA.
//
// In IPA, after multiple folding rounds with challenges
// `ch_1, ch_2, ..., ch_t`, each original element gets multiplied
// by a product of selected challenges depending on its binary index.
// This function precomputes those products.
//
// Example:
// - For challenges [c1, c2], output coefficients are:
//   [1, c2, c1, c1*c2]
pub fn compute_folding_coeffs<C: Curve>(chs: &[C::ScalarField]) -> Vec<C::ScalarField>
where
    <C as Curve>::ScalarField: Arithmetic,
{
    // if no challenges, then panic!
    if chs.is_empty() {
        panic!("empty challenges");
    }
    
    let mut out = vec![C::ScalarField::one()];

    for &ch in chs {
        let n = out.len();
        let mut next = vec![C::ScalarField::zero(); n * 2];

        for j in 0..n {
            let v = out[j];
            next[2 * j] = v;
            next[2 * j + 1] = ch * v;
        }

        out = next;
    }

    out
}