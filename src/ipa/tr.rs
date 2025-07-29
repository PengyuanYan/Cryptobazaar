use icicle_core::curve::{Curve,Affine};
use crate::transcript::TranscriptOracle;
use super::structs::Instance;
use ark_serialize::{CanonicalSerialize, Compress};
use icicle_core::traits::FieldImpl;

pub struct Transcript<const N: usize, const LOG_N: usize, C: Curve> {
    tr: TranscriptOracle<C::ScalarField>,
}

impl<const N: usize, const LOG_N: usize, C: Curve> Transcript<N, LOG_N, C> {
    pub(crate) fn new(init_label: &'static [u8]) -> Self {
        Self {
            tr: TranscriptOracle::new(init_label),
        }
    }

    pub(crate) fn send_instance(&mut self, instance: &Instance<N, C>) {
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();
        self.tr.send_message(b"ipa-instance", &data);
    }

    pub(crate) fn send_l_r(&mut self, l: &Affine::<C>, r: &Affine::<C>) {
        let mut data = Vec::new();
        
        let mut l_x = l.x.to_bytes_le();
        data.append(&mut l_x);
        let mut l_y = l.y.to_bytes_le();
        data.append(&mut l_y);
        self.tr.send_message(b"l", &data);

        let mut r_x = r.x.to_bytes_le();
        data.append(&mut r_x);
        let mut r_y = r.y.to_bytes_le();
        data.append(&mut r_y);
        self.tr.send_message(b"r", &data);
    }

    pub(crate) fn get_alpha_i(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"chi")
    }

    pub(crate) fn get_r(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"r")
    }
}