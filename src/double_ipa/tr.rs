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
        self.tr.send_message(b"double-ipa-instance", &data);
    }

    pub(crate) fn send_ls_rs(
        &mut self,
        l_1: &Affine::<C>,
        r_1: &Affine::<C>,
        l_2: &Affine::<C>,
        r_2: &Affine::<C>,
    ) {
        let mut data = Vec::new();
        
        let mut l_1_x = l_1.x.to_bytes_le();
        data.append(&mut l_1_x);
        let mut l_1_y = l_1.y.to_bytes_le();
        data.append(&mut l_1_y);
        self.tr.send_message(b"l_1", &data);

        let mut r_1_x = r_1.x.to_bytes_le();
        data.append(&mut r_1_x);
        let mut r_1_y = r_1.y.to_bytes_le();
        data.append(&mut r_1_y);
        self.tr.send_message(b"r_1", &data);

        let mut l_2_x = l_2.x.to_bytes_le();
        data.append(&mut l_2_x);
        let mut l_2_y = l_2.y.to_bytes_le();
        data.append(&mut l_2_y);
        self.tr.send_message(b"l_2", &data);

        let mut r_2_x = r_2.x.to_bytes_le();
        data.append(&mut r_2_x);
        let mut r_2_y = r_2.y.to_bytes_le();
        data.append(&mut r_2_y);
        self.tr.send_message(b"r_2", &data);
    }

    pub(crate) fn get_alpha_i(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"chi")
    }

    pub(crate) fn get_r(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"r")
    }
}