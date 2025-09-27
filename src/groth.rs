//! # Groth16 Zero-Knowledge Proof System
//!
//! ## Overview
//!
//! High-level interface for the Groth16 zero-knowledge proof system using arkworks
//! elliptic curve implementations. This module provides a bridge between the generic
//! BasicPairing trait and concrete arkworks pairing implementations.
//!
//! ## 🚨⚠️ **IMPORTANT** ⚠️🚨
//!
//! **THIS IS A LEARNING IMPLEMENTATION - NOT FOR PRODUCTION USE!**
//!
//! While this module uses real cryptographic pairings (BN254, BLS12-381), it lacks
//! critical security features and optimizations.
//!
//! For production systems, use established libraries like arkworks-groth16 or bellman
//! that have undergone security audits and performance optimization.
//!
//! This implementation supports any arkworks `Pairing` implementation.
//!
//! ## Usage Pattern
//!
//! ```rust
//! // 1. Generate QAP from constraints
//! let qap = constraint_system_to_qap(constraints);
//! 
//! // 2. Trusted setup (once per circuit)
//! let setup = trusted_setup::<Bn254>(&qap);
//! 
//! // 3. Generate proof (per instance)
//! let proof = prove::<Bn254>(&qap, &witness, &setup);
//! 
//! // 4. Verify proof (public verification)
//! let valid = verify::<Bn254>(&qap, &public_inputs, &proof, &setup);
//! ```

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::CurveGroup;

use crate::linear_proof::{setup_linear, prove_linear, verify_linear, Groth16SetupParameters, GrothProof};
use crate::basic_pairing::{BasicPairing, BasicPairingGroup};
use crate::quadratic_arithmetic_programs::QAP;

/// Bridge implementation: BasicPairingGroup for arkworks CurveGroup
///
/// Provides a generic implementation of BasicPairingGroup for any arkworks
/// CurveGroup type, enabling compatibility with our abstract pairing interface.
///
/// # Type Parameters
///
/// * `G` - Any type implementing arkworks CurveGroup trait
impl<G> BasicPairingGroup<G::ScalarField> for G
where
    G: CurveGroup,
{
    fn generator() -> Self {
        G::generator()
    }
}

/// Bridge implementation: BasicPairing for arkworks Pairing types
///
/// Provides a generic implementation of BasicPairing for any arkworks pairing
/// type, enabling our abstract algorithms to work with concrete elliptic curve
/// pairings like BN254 and BLS12-381.
///
/// # Type Parameters
///
/// * `P` - Any type implementing arkworks Pairing trait
///
/// # Associated Types
///
/// * `ScalarField` - Scalar field of the elliptic curve
/// * `G1` - First group (typically the base curve)
/// * `G2` - Second group (typically the twisted curve)  
/// * `TargetGroup` - Target group of the pairing operation
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

/// Perform Groth16 trusted setup ceremony for a given QAP
///
/// Generates the proving and verification keys required for the Groth16 protocol
/// using cryptographically secure random toxic waste parameters. This is the
/// initialization phase that must be performed once per constraint system.
///
/// # Type Parameters
///
/// * `P` - Pairing type (e.g., Bn254, Bls12_381)
///
/// # Arguments
///
/// * `qap` - Quadratic arithmetic program representing the constraint system
///
/// # Returns
///
/// Complete setup parameters containing:
/// - sigma1: Proof/verification material to be used with the first group G1
/// - sigma2: Proof/verification material to be used with the second group G2
///
/// # Example
///
/// ```rust
/// use ark_bn254::Bn254;
/// let qap = create_qap_from_constraints(constraints);
/// let setup = trusted_setup::<Bn254>(&qap);
/// // setup can now be used for proving and verification
/// ```
pub fn trusted_setup<P: Pairing>(qap: &QAP<P::ScalarField>) -> Groth16SetupParameters<P::G1, P::G2> {
    setup_linear::<P>(qap)
}

/// Generate a Groth16 zero-knowledge proof
///
/// Creates a succinct zero-knowledge proof that demonstrates knowledge of a
/// satisfying witness for the given QAP without revealing the witness itself.
///
/// # Type Parameters
///
/// * `P` - Pairing type (e.g., Bn254, Bls12_381)
///
/// # Arguments
///
/// * `qap` - Quadratic arithmetic program to prove satisfaction for
/// * `witness` - Complete witness vector (public inputs + private variables)
/// * `setup` - Trusted setup parameters for the circuit
///
/// # Returns
///
/// Three-element proof (A, B, C) that can be publicly verified.
///
/// # Example
///
/// ```rust
/// let witness = vec![public_input, private_var1, private_var2];
/// let proof = prove::<Bn254>(&qap, &witness, &setup);
/// ```
pub fn prove<P: Pairing>(
    qap: &QAP<P::ScalarField>,
    witness: &[P::ScalarField], 
    setup: &Groth16SetupParameters<P::G1, P::G2>
) -> GrothProof<P::G1, P::G2> {
    prove_linear::<P>(qap, witness, setup)
}

/// Verify a Groth16 zero-knowledge proof
///
/// Checks whether a proof is valid for given public inputs using the Groth16
/// verification equation. Verification is fast and requires only public information
/// (including public information from the trusted_setup).
///
/// # Type Parameters
///
/// * `P` - Pairing type (e.g., Bn254, Bls12_381)
///
/// # Arguments
///
/// * `qap` - Quadratic arithmetic program the proof claims to satisfy
/// * `public_inputs` - Public portion of the witness (known to verifier)
/// * `proof` - Three-element proof (A, B, C) to verify
/// * `setup` - Trusted setup parameters (verification key portion)
///
/// # Returns
///
/// * `true` - Proof is valid for the given public inputs
/// * `false` - Proof is invalid or public inputs are incorrect
///
/// # Example
///
/// ```rust
/// let valid = verify::<Bn254>(&qap, &public_inputs, &proof, &setup);
/// if valid {
///     println!("Proof verified successfully!");
/// }
/// ```
pub fn verify<P: Pairing>(
    qap: &QAP<P::ScalarField>,
    public_inputs: &[P::ScalarField],
    proof: &GrothProof<P::G1, P::G2>,
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
