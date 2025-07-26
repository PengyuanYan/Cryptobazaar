use icicle_core::traits::FieldImpl;
use merlin::Transcript as Tr;
use std::marker::PhantomData;
use std::sync::Mutex;

/// Generic structure that serves same interface to all subprotocols
/// Makes it thread safe
pub struct TranscriptOracle<F: FieldImpl> {
    tr: Mutex<Tr>,
    challenge_buffer: Mutex<Vec<u8>>,
    _f: PhantomData<F>,
}

impl<F: FieldImpl> TranscriptOracle<F> {
    pub fn new(init_label: &'static [u8]) -> Self {
        let tr = Tr::new(init_label);
        let challenge_size = F::zero().to_bytes_le().len();
        let challenge_buffer = vec![0u8; challenge_size];
        Self {
            tr: Mutex::new(tr),
            challenge_buffer: Mutex::new(challenge_buffer),
            _f: PhantomData,
        }
    }

    pub fn send_message(&mut self, label: &'static [u8], data: &[u8]) {
        let mut tr = self.tr.lock().unwrap();
        tr.append_message(label, data);
    }

    pub fn squeeze_challenge(&mut self, label: &'static [u8]) -> F {
        let mut tr = self.tr.lock().unwrap();
        let mut ch_buffer = self.challenge_buffer.lock().unwrap();

        tr.challenge_bytes(label, &mut ch_buffer);
        
        // only use the half of the bytes to make sure the from_bytes_le is not overwhelmed
        let half = ch_buffer.len() / 2;
        let first_half: &[u8] = &ch_buffer[..half];
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
    fn test_from_bytes_le() {
        let transcript = Tr::new(b"acc-transcript");
        let challenge_size = Bn254ScalarField::zero().to_bytes_le().len();
        let challenge_buffer = vec![0u8; challenge_size];

        let tr = Mutex::new(transcript);
        let challenge_buffer = Mutex::new(challenge_buffer);

        let mut ch_buffer = challenge_buffer.lock().unwrap();
        tr.lock().unwrap().challenge_bytes(b"acc-transcript", &mut ch_buffer);

        let beta = Bn254ScalarField::from_bytes_le(&ch_buffer);

        let beta_1 = beta + <Bn254CurveCfg as Curve>::ScalarField::zero();
        let beta_2 = beta_1 + <Bn254CurveCfg as Curve>::ScalarField::zero();
        //println!("{}", beta_1 - beta_2);

        let test_coeffs = [<Bn254CurveCfg as Curve>::ScalarField::zero() - beta, <Bn254CurveCfg as Curve>::ScalarField::one()];
        let test_poly = Bn254DensePolynomial::from_coeffs(HostSlice::from_slice(&test_coeffs), 2);
        let x = test_poly.eval(&beta);

        //println!("{}",x);
    }
}