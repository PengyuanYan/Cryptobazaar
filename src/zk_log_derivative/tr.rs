use icicle_core::curve::{Curve,Affine};
use ark_serialize::{CanonicalSerialize, Compress};
use icicle_core::traits::FieldImpl;

use crate::transcript::TranscriptOracle;
use super::structs::{Instance, VerifierIndex};

pub struct Transcript<C: Curve> {
    tr: TranscriptOracle<C::ScalarField>,
}

impl<C: Curve> Transcript<C> {
    pub(crate) fn new(init_label: &'static [u8]) -> Self {
        Self {
            tr: TranscriptOracle::new(init_label),
        }
    }

    pub(crate) fn send_v_index(&mut self, index: &VerifierIndex<C>) {
        let mut data = Vec::new();
        index.serialize_with_mode(&mut data, Compress::No).unwrap();
        self.tr.send_message(b"log-derivative-index", &data);
    }

    pub(crate) fn send_instance(&mut self, instance: &Instance<C>) {
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();
        self.tr.send_message(b"log-derivative-instance", &data);
    }

    pub(crate) fn send_blinders_sum(&mut self, bs: &C::ScalarField) {
        let mut data = Vec::new();
        let mut bs_bytes = bs.to_bytes_le();
        data.append(&mut bs_bytes);
        self.tr.send_message(b"bs", &data);
    }

    pub(crate) fn send_b_and_q(&mut self, b: &Affine::<C>, q: &Affine::<C>) {
        let mut data = Vec::new();

        let mut b_x = b.x.to_bytes_le();
        data.append(&mut b_x);
        let mut b_y = b.y.to_bytes_le();
        data.append(&mut b_y);
        self.tr.send_message(b"b", &data);

        let mut q_x = q.x.to_bytes_le();
        data.append(&mut q_x);
        let mut q_y = q.y.to_bytes_le();
        data.append(&mut q_y);
        self.tr.send_message(b"q", &data);
    }

    pub(crate) fn send_openings(
        &mut self,
        f_opening: &C::ScalarField,
        s_opening: &C::ScalarField,
        b_opening: &C::ScalarField,
        q_opening: &C::ScalarField,
    ) {
        let mut data = Vec::new();
        
        let mut f_opening_bytes = f_opening.to_bytes_le();
        data.append(&mut f_opening_bytes);
        self.tr.send_message(b"f_opening", &data);

        let mut s_opening_bytes = s_opening.to_bytes_le();
        data.append(&mut s_opening_bytes);
        self.tr.send_message(b"s_opening", &data);

        let mut b_opening_bytes = b_opening.to_bytes_le();
        data.append(&mut b_opening_bytes);
        self.tr.send_message(b"b_opening", &data);

        let mut q_opening_bytes = q_opening.to_bytes_le();
        data.append(&mut q_opening_bytes);
        self.tr.send_message(b"q_opening", &data);
    }

    pub(crate) fn get_beta(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"beta")
    }

    pub(crate) fn get_mu(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"mu")
    }

    pub(crate) fn get_separation_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"separation-challenge")
    }
}