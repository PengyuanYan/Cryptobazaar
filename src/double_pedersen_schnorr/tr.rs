use icicle_core::curve::{Curve,Affine};
use crate::transcript::TranscriptOracle;
use super::structs::Instance;
use ark_serialize::{CanonicalSerialize, Compress};
use icicle_core::traits::FieldImpl;

pub struct Transcript<C: Curve> {
    tr: TranscriptOracle<C::ScalarField>,
}

impl<C: Curve> Transcript<C> {
    pub(crate) fn new(init_label: &'static [u8]) -> Self {
        Self {
            tr: TranscriptOracle::new(init_label),
        }
    }

    pub(crate) fn send_instance(&mut self, instance: &Instance<C>) {
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();
        self.tr
            .send_message(b"double-pedersen-schnorr-instance", &data);
    }

    pub(crate) fn send_blinders(&mut self, rand_1: &Affine::<C>, rand_2: &Affine::<C>) {
        let mut data = Vec::new();
        
        let mut rand_1_x = rand_1.x.to_bytes_le();
        data.append(&mut rand_1_x);
        let mut rand_1_y = rand_1.y.to_bytes_le();
        data.append(&mut rand_1_y);
        self.tr.send_message(b"r", &data);

        let mut rand_2_x = rand_2.x.to_bytes_le();
        data.append(&mut rand_2_x);
        let mut rand_2_y = rand_2.y.to_bytes_le();
        data.append(&mut rand_2_y);
        self.tr.send_message(b"q", &data);
    }

    pub(crate) fn get_c(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"c")
    }
}