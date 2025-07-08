use icicle_core::curve::Curve;

use super::structs::Instance;
use crate::transcript::TranscriptOracle;

use icicle_core::traits::FieldImpl;

use ark_serialize::{CanonicalSerialize, Compress};
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
        self.tr.send_message(b"acc-instance", &data);
    }

    pub(crate) fn send_q(&mut self, q: &C::ScalarField) {
        let mut data = q.to_bytes_le();
        self.tr.send_message(b"q", &data);
    }

    pub(crate) fn get_beta(&mut self) -> C::ScalarField {
        self.tr.squeeze_challenge(b"beta")
    }
}

#[cfg(test)]
mod tr_tests {
    use super::Transcript;
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;

    #[test]
    fn test_tr() {
        let mut tr1: Transcript<Bn254CurveCfg> = Transcript::new(b"acc-transcript");
        let beta_1 = tr1.get_beta();

        let mut tr2 = Transcript::<Bn254CurveCfg>::new(b"acc-transcript");
        let beta_2 = tr2.get_beta();
        assert_eq!(beta_1,beta_2);
    }
}