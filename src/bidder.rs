use icicle_core::curve::{Curve,Affine};
use icicle_core::pairing::Pairing;
use icicle_core::traits::FieldImpl;
use icicle_core::polynomials::UnivariatePolynomial;
use icicle_core::ntt::{NTT};
use icicle_core::traits::Arithmetic;
use std::ops::{Add, Mul};
use rand::{RngCore, SeedableRng};

use crate::{
    bid_encoder::BidEncoder,
    gates::{
        structs::{
            Proof as GatesProof, ProverIndex as GProverIndex, VerifierIndex as GVerifierIndex,
            Witness as GatesWitness,
        },
        GatesArgument,
    },
    kzg::PK as KzgPk,
};

pub struct Bidder<const P: usize, const N: usize, C1, C2, F, U>
where
    C1: Curve,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = C1::ScalarField>,
{
    pk: KzgPk<C1, C2, F>,
    gp_index: GProverIndex<C1, U>,
    gv_index: GVerifierIndex<C1>,
    bid_encoder: Option<BidEncoder<P, N, C1>>,
}

impl<const P: usize, const N: usize, C1, C2, F, U> Bidder<P, N, C1, C2, F, U>
where
    C1: Curve + icicle_core::msm::MSM<C1>,
    C2: Curve,
    F: FieldImpl,
    C1: Pairing<C1, C2, F>,
    U: UnivariatePolynomial<Field = C1::ScalarField>,
    for<'a> &'a U: Add<&'a U, Output = U>
{
    pub fn new(pk: KzgPk<C1, C2, F>) -> Self 
    where
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic
    {
        let gp_index = GatesArgument::<N, P, C1, C2, F>::prover_index();
        let gv_index = GatesArgument::<N, P, C1, C2, F>::verifier_index::<U>(&pk);
        Self {
            pk,
            gp_index,
            gv_index,
            bid_encoder: None,
        }
    }
    
    pub fn encode<R: RngCore + SeedableRng>(&mut self, bid: usize, seed: R::Seed) {
        self.bid_encoder = Some(BidEncoder::encode::<R>(bid, seed));
    }

    pub fn construct_bid_well_formation_proof<R: RngCore + SeedableRng>(
        &self,
        seed: R::Seed,
    ) -> GatesProof<C1> 
    where
        <<C1 as Curve>::ScalarField as FieldImpl>::Config: NTT<<C1 as Curve>::ScalarField, <C1 as Curve>::ScalarField>,
        <C1 as Curve>::ScalarField: Arithmetic
    {
        let bid_encoder = self.bid_encoder.as_ref().unwrap();
        let witness: GatesWitness<C1, U> = bid_encoder.to_gate_witness::<R, U>(seed);
        GatesArgument::<N, P, C1, C2, F>::prove(&witness, &self.gv_index, &self.gp_index, &self.pk)
    }

    pub fn first_round(&self) -> Vec<Affine::<C1>> {
        let bid_encoder = self.bid_encoder.as_ref().unwrap();
        bid_encoder.to_first_av_round()
    }

    pub fn second_round(&self, basis: &[Affine::<C1>]) -> Vec<Affine::<C1>>
    where 
        <C1 as Curve>::ScalarField: Add<Output = <C1 as Curve>::ScalarField> + Mul<Output = <C1 as Curve>::ScalarField>
    {
        let bid_encoder = self.bid_encoder.as_ref().unwrap();
        bid_encoder.to_second_av_round(&basis)
    }
}