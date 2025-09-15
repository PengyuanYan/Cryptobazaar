// This is the Oracle of Fiat-Shamir Transformaion for univariate sumcheck.
use icicle_core::curve::{Curve,Affine};
use crate::transcript::TranscriptOracle;
use super::structs::Instance;
use ark_serialize::{CanonicalSerialize, Compress};
use icicle_core::traits::FieldImpl;

pub struct Transcript<C: Curve> {
    tr: TranscriptOracle<C::ScalarField>,
}

impl<C: Curve> Transcript<C> {
    pub(crate) fn new_transcript(init_label: &'static [u8]) -> Self {
        Self {
            tr: TranscriptOracle::new_transcript(init_label),
        }
    }

    pub(crate) fn send_instance(&mut self, instance: &Instance<C>) {
        let mut data = Vec::new();
        instance.serialize_with_mode(&mut data, Compress::No).unwrap();
        self.tr.send_message(b"univariate-sumcheck-instance", &data);
    }

    pub(crate) fn send_r_and_q(&mut self, r_cm: &Affine::<C>, q_cm: &Affine::<C>) {
        let mut data = Vec::new();

        let mut r_cm_x = r_cm.x.to_bytes_le();
        data.append(&mut r_cm_x);
        let mut r_cm_y = r_cm.y.to_bytes_le();
        data.append(&mut r_cm_y);
        self.tr.send_message(b"r", &data);
        
        let mut q_cm_x = q_cm.x.to_bytes_le();
        data.append(&mut q_cm_x);
        let mut q_cm_y = q_cm.y.to_bytes_le();
        data.append(&mut q_cm_y);
        self.tr.send_message(b"q", &data);
    }

    pub(crate) fn send_openings(
        &mut self,
        a_eval: &C::ScalarField,
        b_eval: &C::ScalarField,
        r_eval: &C::ScalarField,
        q_eval: &C::ScalarField,
    ) {
        let mut scalars_data = Vec::new();
        

        let mut a_eval_bytes = a_eval.to_bytes_le();
        scalars_data.append(&mut a_eval_bytes);
        self.tr.send_message(b"a_eval", &scalars_data);

        let mut b_eval_bytes = b_eval.to_bytes_le();
        scalars_data.append(&mut b_eval_bytes);
        self.tr.send_message(b"b_eval", &scalars_data);

        let mut r_eval_bytes = r_eval.to_bytes_le();
        scalars_data.append(&mut r_eval_bytes);
        self.tr.send_message(b"r_eval", &scalars_data);

        let mut q_eval_bytes = q_eval.to_bytes_le();
        scalars_data.append(&mut q_eval_bytes);
        self.tr.send_message(b"q_eval", &scalars_data);
    }

    pub(crate) fn get_opening_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"opening_challenge")
    }

    pub(crate) fn get_separation_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"separation_challenge")
    }
}