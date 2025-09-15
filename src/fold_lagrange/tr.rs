// This is the Oracle of Fiat-Shamir Transformaion for fold_lagrange.
use icicle_core::curve::{Curve,Affine};
use icicle_core::traits::FieldImpl;

use crate::transcript::TranscriptOracle;

pub struct Transcript<C: Curve> {
    tr: TranscriptOracle<C::ScalarField>,
}

impl<C: Curve> Transcript<C> {
    pub(crate) fn new_transcript(init_label: &'static [u8]) -> Self {
        Self {
            tr: TranscriptOracle::new_transcript(init_label),
        }
    }

    pub(crate) fn send_p(&mut self, p: &Affine::<C>) {
        let mut data = Vec::new();
        
        let mut p_x = p.x.to_bytes_le();
        data.append(&mut p_x);
        let mut p_y = p.y.to_bytes_le();
        data.append(&mut p_y);
        self.tr.send_message(b"p", &data);
    }

    pub(crate) fn get_beta(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"beta")
    }
}