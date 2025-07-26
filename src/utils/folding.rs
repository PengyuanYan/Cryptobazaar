use icicle_core::curve::{Curve, Affine};
use icicle_core::traits::FieldImpl;
use icicle_core::traits::Arithmetic;

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

        let left = &x[..x.len() / 2];
        let right = &x[x.len() / 2..];

        let mut folded = Vec::with_capacity(left.len());

        for i in 0..left.len() {
            let l = left[i];
            let r = right[i];

            folded.push(l + r * ch);
        }

        Ok(folded)
    }
}

pub struct AffFold<C: Curve>(C);

impl<C: Curve> Fold for AffFold<C> {
    type Challenge = C::ScalarField;
    type FoldType = Affine::<C>;
    fn fold_vec(
        x: &[Self::FoldType],
        ch: Self::Challenge,
    ) -> Result<Vec<Self::FoldType>, FoldError> {
        Self::check(x)?;

        let left = &x[..x.len() / 2];
        let right = &x[x.len() / 2..];

        let affine_result: Vec<Affine::<C>> = {
            let mut tmp = Vec::with_capacity(left.len());
            for i in 0..left.len() {
                let result_projective = left[i].to_projective() + C::mul_scalar(right[i].to_projective(), ch);

                let mut result_affine = Affine::<C>::zero();
                C::to_affine(&result_projective, &mut result_affine);

                tmp.push(result_affine);
            }
            tmp
        };

        Ok(affine_result)
    }
}

pub fn compute_folding_coeffs<C: Curve>(chs: &[C::ScalarField]) -> Vec<C::ScalarField>
where
    <C as Curve>::ScalarField: Arithmetic,
{
    let l = chs.len();

    let mut c: Vec<Vec<C::ScalarField>> = Vec::with_capacity(l);
    for i in 0..l {
        c.push(vec![C::ScalarField::one(); 1 << (i + 1)])
    }

    // first array is equal to [1][ch_0]
    c[0][1] = chs[0];

    for i in 0..(l - 1) {
        for j in 0..(1 << (i + 1)) {
            c[i + 1][2 * j] = c[i][j];
            c[i + 1][2 * j + 1] = chs[i + 1] * c[i][j];
        }
    }

    c[l - 1].clone()
}