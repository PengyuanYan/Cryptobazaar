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
        self.tr.send_message(b"pedersen-schnorr-instance", &data);
    }

    pub(crate) fn send_blinder(&mut self, blinder: &Affine::<C>) {
        let mut data = Vec::new();

        let mut blinder_x = blinder.x.to_bytes_le();
        data.append(&mut blinder_x);
        let mut blinder_y = blinder.y.to_bytes_le();
        data.append(&mut blinder_y);
        self.tr.send_message(b"blinder", &data);
    }

    pub(crate) fn get_c(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"c")
    }
}