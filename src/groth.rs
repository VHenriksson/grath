use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::CurveGroup;

use crate::linear_proof::{setup_linear, prove_linear, verify_linear, Groth16SetupParameters, NILPProof};
use crate::pairing_traits::{BasicPairing, BasicPairingGroup};
use crate::quadratic_arithmetic_programs::QAP;

// Implement BasicPairingGroup for any CurveGroup that meets our requirements
impl<G> BasicPairingGroup<G::ScalarField> for G
where
    G: CurveGroup,
{
    fn generator() -> Self {
        G::generator()
    }
}

// Implement BasicPairing for any arkworks Pairing type
impl<P> BasicPairing for P
where
    P: Pairing,
{
    type ScalarField = P::ScalarField;
    type G1 = P::G1;
    type G2 = P::G2;
    type TargetGroup = PairingOutput<P>;

    fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetGroup {
        P::pairing(g1_elem, g2_elem)
    }
}

/// Groth16 trusted setup - generates the proving and verification keys
pub fn trusted_setup<P: Pairing>(qap: &QAP<P::ScalarField>) -> Groth16SetupParameters<P::G1, P::G2> {
    setup_linear::<P>(qap)
}

/// Groth16 prove function - generates a zero-knowledge proof
pub fn prove<P: Pairing>(
    qap: &QAP<P::ScalarField>,
    witness: &[P::ScalarField], 
    setup: &Groth16SetupParameters<P::G1, P::G2>
) -> NILPProof<P::G1, P::G2> {
    prove_linear::<P>(qap, witness, setup)
}

/// Groth16 verify function - verifies a zero-knowledge proof
pub fn verify<P: Pairing>(
    qap: &QAP<P::ScalarField>,
    public_inputs: &[P::ScalarField],
    proof: &NILPProof<P::G1, P::G2>,
    setup: &Groth16SetupParameters<P::G1, P::G2>
) -> bool {
    verify_linear::<P>(qap, public_inputs, proof, setup)
}

#[cfg(test)]
mod tests {

    use crate::test_functions::{complex_qap, simple_qap};

    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_bn254::{Bn254, Fr};
    use ark_ec::bls12::Bls12;
    use ark_ff::{Field, Zero, One};
    use ark_poly::univariate::DensePolynomial;


    #[test]
    fn test_groth16_with_bn254() {
        let qap = simple_qap::qap();
        let setup = trusted_setup::<Bn254>(&qap);
        let proof = prove::<Bn254>(&qap, &simple_qap::witness(), &setup);
        assert!(verify::<Bn254>(&qap, &simple_qap::public_input(), &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify::<Bn254>(&qap, &simple_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
    }


    #[test]
    fn test_complex_qap() {
        let qap = complex_qap::qap::<Fr>();
        let setup = trusted_setup::<Bn254>(&qap);
        let proof = prove::<Bn254>(&qap, &complex_qap::witness(), &setup);
        assert!(verify::<Bn254>(&qap, &complex_qap::public_input(), &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify::<Bn254>(&qap, &complex_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
    }

    #[test]
    fn test_complex_qap_different_witness() {
        let qap = complex_qap::qap::<Fr>();
        let setup = trusted_setup::<Bn254>(&qap);
        let proof = prove::<Bn254>(&qap, &complex_qap::witness_alternative(), &setup);
        assert!(verify::<Bn254>(&qap, &complex_qap::public_input_alternative(), &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify::<Bn254>(&qap, &complex_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
    }

    #[test]
    fn test_complex_qap_with_different_curve() {
        use ark_bls12_381::Fr;
        let qap = complex_qap::qap::<Fr>();
        let setup = trusted_setup::<Bls12_381>(&qap);
        let proof = prove::<Bls12_381>(&qap, &complex_qap::witness(), &setup);
        assert!(verify::<Bls12_381>(&qap, &complex_qap::public_input(), &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify::<Bls12_381>(&qap, &complex_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
    }

}
