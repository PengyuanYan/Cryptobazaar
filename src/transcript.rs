use icicle_core::traits::FieldImpl;
use merlin::Transcript as Tr;
use std::marker::PhantomData;
use std::sync::Mutex;

// This module defines a generic, thread-safe wrapper around the Merlin transcript
// The code is directly used from the original Cryptobazaar.
pub struct TranscriptOracle<F: FieldImpl> {
    transcript: Mutex<Tr>,
    challenge_buffer: Mutex<Vec<u8>>,
    _f: PhantomData<F>,
}

impl<F: FieldImpl> TranscriptOracle<F> {
    /// Create a new transcript with an initialization label.
    pub fn new_transcript(init_label: &'static [u8]) -> Self {
        let transcript = Tr::new(init_label);
        let challenge_size = F::zero().to_bytes_le().len();
        let challenge_buffer = vec![0u8; challenge_size];
        Self {
            transcript: Mutex::new(transcript),
            challenge_buffer: Mutex::new(challenge_buffer),
            _f: PhantomData,
        }
    }
    
    // Append a message to the transcript under a given label.
    // Used to commit public data into the Fiat–Shamir transcript.
    pub fn send_message(&mut self, label: &'static [u8], data: &[u8]) {
        let mut transcript = self.transcript.lock().unwrap();
        transcript.append_message(label, data);
    }
    
    // Sample a challenge from the transcript, labeled by `label`.
    // The challenge is deterministically derived from the transcript state.
    pub fn squeeze_challenge(&mut self, label: &'static [u8]) -> F {
        let mut transcript = self.transcript.lock().unwrap();
        let mut challenge_buffer = self.challenge_buffer.lock().unwrap();

        transcript.challenge_bytes(label, &mut challenge_buffer);
        
        // only use the half of the bytes to make sure the from_bytes_le is not overwhelmed
        let half = challenge_buffer.len() / 2;
        let first_half: &[u8] = &challenge_buffer[..half];
        F::from_bytes_le(&first_half)
    }
}

#[cfg(test)]
mod tr_tests {
    use merlin::Transcript as Tr;
    use icicle_bn254::curve::CurveCfg as Bn254CurveCfg;
    use icicle_bn254::curve::ScalarField as Bn254ScalarField;
    use icicle_core::traits::FieldImpl;
    use icicle_core::curve::Curve;
    use std::sync::Mutex;
    use icicle_bn254::polynomials::DensePolynomial as Bn254DensePolynomial;
    use icicle_core::polynomials::UnivariatePolynomial;
    use icicle_runtime::memory::HostSlice;

    #[test]
    #[ignore]
    fn test_from_bytes_le() {
        let transcript = Tr::new(b"acc-transcript");
        let challenge_size = Bn254ScalarField::zero().to_bytes_le().len();
        let challenge_buffer = vec![0u8; challenge_size];

        let tr = Mutex::new(transcript);
        let challenge_buffer = Mutex::new(challenge_buffer);

        let mut ch_buffer = challenge_buffer.lock().unwrap();
        tr.lock().unwrap().challenge_bytes(b"acc-transcript", &mut ch_buffer);

        let beta = Bn254ScalarField::from_bytes_le(&ch_buffer);

        let _beta_1 = beta + <Bn254CurveCfg as Curve>::ScalarField::zero();
        let _beta_2 = _beta_1 + <Bn254CurveCfg as Curve>::ScalarField::zero();
        //println!("{}", beta_1 - beta_2);

        let test_coeffs = [<Bn254CurveCfg as Curve>::ScalarField::zero() - beta, <Bn254CurveCfg as Curve>::ScalarField::one()];
        let test_poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&test_coeffs), 2);
        let _x = test_poly.eval(&beta);

        //println!("{}",x);
    }
}