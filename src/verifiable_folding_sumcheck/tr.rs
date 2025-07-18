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
        self.tr.send_message(b"verifiable-folding-instance", &data);
    }

    pub(crate) fn send_blinders(&mut self, s: &Affine::<C>, blinder: &Affine::<C>) {
        let mut data = Vec::new();
        
        let mut s_x = s.x.to_bytes_le();
        data.append(&mut s_x);
        let mut s_y = s.y.to_bytes_le();
        data.append(&mut s_y);
        self.tr.send_message(b"s", &data);
        
        let mut blinder_x = blinder.x.to_bytes_le();
        data.append(&mut blinder_x);
        let mut blinder_y = blinder.y.to_bytes_le();
        data.append(&mut blinder_y);
        self.tr.send_message(b"blinder", &data);
    }

    pub(crate) fn second_round(
        &mut self,
        z_1: &C::ScalarField,
        z_2: &C::ScalarField,
        r_cm: &Affine::<C>,
        r_degree_cm: &Affine::<C>,
        q_cm: &Affine::<C>,
    ) {
        let mut points_data = Vec::new();
        let mut scalars_data = Vec::new();

        let mut r_cm_x = r_cm.x.to_bytes_le();
        points_data.append(&mut r_cm_x);
        let mut r_cm_y = r_cm.y.to_bytes_le();
        points_data.append(&mut r_cm_y);
        self.tr.send_message(b"r", &points_data);

        let mut r_degree_cm_x = r_degree_cm.x.to_bytes_le();
        points_data.append(&mut r_degree_cm_x);
        let mut r_degree_cm_y = r_degree_cm.y.to_bytes_le();
        points_data.append(&mut r_degree_cm_y);
        self.tr.send_message(b"r_degree", &points_data);

        let mut q_cm_x = q_cm.x.to_bytes_le();
        points_data.append(&mut q_cm_x);
        let mut q_cm_y = q_cm.y.to_bytes_le();
        points_data.append(&mut q_cm_y);
        self.tr.send_message(b"q", &points_data);

        let mut z_1_bytes = z_1.to_bytes_le();
        scalars_data.append(&mut z_1_bytes);
        self.tr.send_message(b"z_1", &scalars_data);

        let mut z_2_bytes = z_2.to_bytes_le();
        scalars_data.append(&mut z_2_bytes);
        self.tr.send_message(b"z_2", &scalars_data);
    }

    pub(crate) fn send_openings(
        &mut self,
        a_opening: &C::ScalarField,
        blinder_opening: &C::ScalarField,
        r_opening: &C::ScalarField,
        q_opening: &C::ScalarField,
    ) {
        let mut data = Vec::new();
        
        let mut a_opening_bytes = a_opening.to_bytes_le();
        data.append(&mut a_opening_bytes);
        self.tr.send_message(b"a_opening", &data);

        let mut blinder_opening_bytes = blinder_opening.to_bytes_le();
        data.append(&mut blinder_opening_bytes);
        self.tr.send_message(b"blinder_opening", &data);

        let mut r_opening_bytes = r_opening.to_bytes_le();
        data.append(&mut r_opening_bytes);
        self.tr.send_message(b"r_opening", &data);

        let mut q_opening_bytes = q_opening.to_bytes_le();
        data.append(&mut q_opening_bytes);
        self.tr.send_message(b"q_opening", &data);
    }

    pub(crate) fn get_c(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"c")
    }

    pub(crate) fn get_opening_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"opening_challenge")
    }

    pub(crate) fn get_separation_challenge(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"separation_challenge")
    }
}