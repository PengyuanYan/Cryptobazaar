// This is the Oracle of Fiat-Shamir Transformaion for Acc
use icicle_core::curve::{Curve,Affine};
use crate::transcript::TranscriptOracle;
use super::structs::VerifierIndex;
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

    pub(crate) fn send_index(&mut self, index: &VerifierIndex<C>) {
        let mut data = Vec::new();
        index.serialize_with_mode(&mut data, Compress::No).unwrap();
        self.tr.send_message(b"gates-v-index", &data);
    }

    pub(crate) fn send_oracle_commitments(
        &mut self,
        bid: &Affine::<C>,
        f: &Affine::<C>,
        r: &Affine::<C>,
        r_inv: &Affine::<C>,
        diff: &Affine::<C>,
        g: &Affine::<C>,
    ) {
        let mut data = Vec::new();

        let mut bid_x = bid.x.to_bytes_le();
        data.append(&mut bid_x);
        let mut bid_y = bid.y.to_bytes_le();
        data.append(&mut bid_y);
        self.tr.send_message(b"bid", &data);

        let mut f_x = f.x.to_bytes_le();
        data.append(&mut f_x);
        let mut f_y = f.y.to_bytes_le();
        data.append(&mut f_y);
        self.tr.send_message(b"f", &data);

        let mut r_x = r.x.to_bytes_le();
        data.append(&mut r_x);
        let mut r_y = r.y.to_bytes_le();
        data.append(&mut r_y);
        self.tr.send_message(b"r", &data);

        let mut r_inv_x = r_inv.x.to_bytes_le();
        data.append(&mut r_inv_x);
        let mut r_inv_y = r_inv.y.to_bytes_le();
        data.append(&mut r_inv_y);
        self.tr.send_message(b"r_inv", &data);

        let mut diff_x = diff.x.to_bytes_le();
        data.append(&mut diff_x);
        let mut diff_y = diff.y.to_bytes_le();
        data.append(&mut diff_y);
        self.tr.send_message(b"diff", &data);

        let mut g_x = g.x.to_bytes_le();
        data.append(&mut g_x);
        let mut g_y = g.y.to_bytes_le();
        data.append(&mut g_y);
        self.tr.send_message(b"g", &data);
    }

    pub(crate) fn send_q_chunks(&mut self, q_chunk_0: &Affine::<C>, q_chunk_1: &Affine::<C>) {
        let mut data = Vec::new();

        let q_chunk_0_x = q_chunk_0.x.to_bytes_le();
        data.append(&mut q_chunk_0_x.clone());
        let q_chunk_0_y = q_chunk_0.y.to_bytes_le();
        data.append(&mut q_chunk_0_y.clone());
        self.tr.send_message(b"q_chunk_0", &data);

        let q_chunk_1_x = q_chunk_1.x.to_bytes_le();
        data.append(&mut q_chunk_1_x.clone());
        let q_chunk_1_y = q_chunk_1.y.to_bytes_le();
        data.append(&mut q_chunk_1_y.clone());
        self.tr.send_message(b"q_chunk_1", &data);
    }

    pub(crate) fn send_oracle_openings(
        &mut self,
        q_price: &C::ScalarField,
        bid: &C::ScalarField,
        bid_shift: &C::ScalarField,
        f: &C::ScalarField,
        r: &C::ScalarField,
        r_inv: &C::ScalarField,
        diff: &C::ScalarField,
        g: &C::ScalarField,
        q_chunk_0_opening: &C::ScalarField,
        q_chunk_1_opening: &C::ScalarField,
    ) {
        let mut data = Vec::new();
        
        let mut q_price_bytes = q_price.to_bytes_le();
        data.append(&mut q_price_bytes);
        self.tr.send_message(b"q_price", &data);

        let mut bid_eval_bytes = bid.to_bytes_le();
        data.append(&mut bid_eval_bytes);
        self.tr.send_message(b"bid_eval", &data);

        let mut bid_shift_bytes = bid_shift.to_bytes_le();
        data.append(&mut bid_shift_bytes);
        self.tr.send_message(b"bid_shift", &data);

        let mut f_eval_bytes = f.to_bytes_le();
        data.append(&mut f_eval_bytes);
        self.tr.send_message(b"f_eval", &data);

        let mut r_eval_bytes = r.to_bytes_le();
        data.append(&mut r_eval_bytes);
        self.tr.send_message(b"r_eval", &data);

        let mut r_inv_eval_bytes = r_inv.to_bytes_le();
        data.append(&mut r_inv_eval_bytes);
        self.tr.send_message(b"r_inv_eval", &data);

        let mut diff_eval_bytes = diff.to_bytes_le();
        data.append(&mut diff_eval_bytes);
        self.tr.send_message(b"diff_eval", &data);

        let mut g_eval_bytes = g.to_bytes_le();
        data.append(&mut g_eval_bytes);
        self.tr.send_message(b"g_eval", &data);

        let mut q_chunk_0_opening_bytes = q_chunk_0_opening.to_bytes_le();
        data.append(&mut q_chunk_0_opening_bytes);
        self.tr.send_message(b"q_chunk_0_opening", &data);

        let mut q_chunk_1_opening_bytes = q_chunk_1_opening.to_bytes_le();
        data.append(&mut q_chunk_1_opening_bytes);
        self.tr.send_message(b"q_chunk_1_opening", &data);
    }

    pub(crate) fn get_quotient_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"quotient_challenge")
    }

    pub(crate) fn get_evaluation_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"evaluation_challenge")
    }

    pub(crate) fn get_separation_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"separation_challenge")
    }
}