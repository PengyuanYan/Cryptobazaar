use icicle_core::curve::Curve;
use icicle_core::traits::FieldImpl;
use std::ops::Mul;

use super::is_pow_2;

#[derive(Debug)]
pub enum FoldError {
    TooShort,
    NotPowTwo,
}

pub fn compute_folding_coeffs<C: Curve>(chs: &[C::ScalarField]) -> Vec<C::ScalarField>
where
    <C as Curve>::ScalarField: Mul<Output = C::ScalarField>
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