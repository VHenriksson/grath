use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::CurveGroup;

use crate::linear_proof::{BasicPairing, BasicPairingGroup, setup_linear, prove_linear, verify_linear, Groth16SetupParameters, NILPProof};
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

mod tests {

    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_bn254::{Bn254, Fr};
    use ark_ec::bls12::Bls12;
    use ark_ff::{Field, Zero, One};
    use ark_poly::univariate::DensePolynomial;

    /// Create a QAP that represents the constraint: x * 1 = 4
    /// Variables: [1, x] where:
    /// - variable 0 is the constant 1 (public input)
    /// - variable 1 is x (witness variable)
    fn create_test_qap<F: Field>() -> QAP<F> {
        use ark_poly::DenseUVPolynomial;
        
        // u polynomial for "x * 1": u_0=0, u_1=1 (selects variable 1 = x)
        let u_poly_0 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // constant 0
        let u_poly_1 = DensePolynomial::from_coefficients_vec(vec![F::one()]);  // constant 1 (selects x)
        
        // v polynomial for "x * 1": v_0=1, v_1=0 (selects variable 0 = 1)  
        let v_poly_0 = DensePolynomial::from_coefficients_vec(vec![F::one()]);  // constant 1 (selects 1)
        let v_poly_1 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // constant 0

        // w polynomial for "= 4": w_0=4, w_1=0 (constant 4)
        let w_poly_0 = DensePolynomial::from_coefficients_vec(vec![F::from(4u32)]); // constant 4
        let w_poly_1 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // constant 0

        // Target polynomial: (x-1) - evaluates to 0 when x=1
        let target_poly = DensePolynomial::from_coefficients_vec(vec![-F::one(), F::one()]); // (x-1)

        QAP {
            num_variables: 2,  
            num_inputs: 1,     // public input: the constant 1
            u_polynomials: vec![u_poly_0, u_poly_1],
            v_polynomials: vec![v_poly_0, v_poly_1],
            w_polynomials: vec![w_poly_0, w_poly_1],
            target_polynomial: target_poly,
        }
    }


    #[test]
    fn test_groth16_with_bn254() {
        let qap = create_test_qap();

        // Run trusted setup
        let setup = trusted_setup::<Bn254>(&qap);
        

        // Create a witness: [1, 5, 25] representing x=5, x^2=25
        let witness = vec![
            Fr::one(),        // variable 0: constant 1 (public input)
            Fr::from(4u32),   // variable 1: x = 4 (witness)
        ];
        
        // Public inputs: [1] (just the constant)
        let public_inputs = vec![Fr::one()];
        
        // Generate proof
        let proof = prove::<Bn254>(&qap, &witness, &setup);
        
        // Verify proof
        let is_valid = verify::<Bn254>(&qap, &public_inputs, &proof, &setup);
        assert!(is_valid, "Proof should be valid");
        
        // Test with wrong public input - should fail
        let wrong_public_inputs = vec![Fr::from(2u32)];
        let is_invalid = verify::<Bn254>(&qap, &wrong_public_inputs, &proof, &setup);
        assert!(!is_invalid, "Proof should be invalid with wrong public input");
    }

    /// Create a complex QAP for testing:
    /// Constraint 1: x * x = x2          (x^2)
    /// Constraint 2: x2 * x = x3         (x^3) 
    /// Constraint 3: a * x3 = ax3        (a*x^3)
    /// Constraint 4: b * x2 = bx2        (b*x^2)
    /// Constraint 5: c * x = cx          (c*x)
    /// Constraint 6: 1 * ax3 = ax3       (copy constraint for ax3)
    /// Constraint 7: 1 * bx2 = bx2       (copy constraint for bx2)  
    /// Constraint 8: 1 * cx = cx         (copy constraint for cx)
    /// Constraint 9: 1 * d = d           (copy constraint for d)
    /// Constraint 10: (ax3+bx2+cx+d) * 1 = result (final sum - this encodes addition as single constraint)
    ///
    /// Variables: [1, x, result, a, b, c, d, x2, x3, ax3, bx2, cx] (12 variables)
    /// Public inputs: constant 1, x, and result (first 3 variables as required)
    fn create_complex_qap<F: Field>() -> QAP<F> {
        use ark_relations::r1cs::ConstraintMatrices;
        use crate::quadratic_arithmetic_programs::constraint_matrices_to_qap;
        
        // Variables: [1, x, result, a, b, c, d, x2, x3, ax3, bx2, cx]
        //           [ 0, 1,      2, 3, 4, 5, 6,  7,  8,   9,  10, 11]
        
        let constraint_matrices = ConstraintMatrices {
            num_instance_variables: 1,  // Just the constant 1
            num_witness_variables: 11,  // x, result, a, b, c, d, x2, x3, ax3, bx2, cx
            num_constraints: 6,
            a: vec![
                vec![(F::one(), 1)],        // Constraint 0: x * x = x2
                vec![(F::one(), 7)],        // Constraint 1: x2 * x = x3
                vec![(F::one(), 3)],        // Constraint 2: a * x3 = ax3
                vec![(F::one(), 4)],        // Constraint 3: b * x2 = bx2
                vec![(F::one(), 5)],        // Constraint 4: c * x = cx
                vec![(F::one(), 9), (F::one(), 10), (F::one(), 11), (F::one(), 6)],  // Constraint 5: (ax3+bx2+cx+d) * 1 = result
            ],
            a_num_non_zero: 8,
            b: vec![
                vec![(F::one(), 1)],        // Constraint 0: x * x = x2
                vec![(F::one(), 1)],        // Constraint 1: x2 * x = x3
                vec![(F::one(), 8)],        // Constraint 2: a * x3 = ax3
                vec![(F::one(), 7)],        // Constraint 3: b * x2 = bx2
                vec![(F::one(), 1)],        // Constraint 4: c * x = cx
                vec![(F::one(), 0)],        // Constraint 5: (ax3+bx2+cx+d) * 1 = result
            ],
            b_num_non_zero: 6,
            c: vec![
                vec![(F::one(), 7)],        // Constraint 0: x * x = x2
                vec![(F::one(), 8)],        // Constraint 1: x2 * x = x3
                vec![(F::one(), 9)],        // Constraint 2: a * x3 = ax3
                vec![(F::one(), 10)],       // Constraint 3: b * x2 = bx2
                vec![(F::one(), 11)],       // Constraint 4: c * x = cx
                vec![(F::one(), 2)],        // Constraint 5: (ax3+bx2+cx+d) * 1 = result
            ],
            c_num_non_zero: 6,
        };
        
        constraint_matrices_to_qap(&constraint_matrices, 3)
    }

    #[test]
    fn test_complex_qap() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_complex_qap::<Fr>();
        

        let setup = trusted_setup::<Bn254>(&qap);

        // Create a satisfying witness for polynomial: result = a*x^3 + b*x^2 + c*x + d
        // Let's use: x=3, a=2, b=1, c=4, d=5
        // Expected: result = 2*27 + 1*9 + 4*3 + 5 = 54 + 9 + 12 + 5 = 80
        
        let x = Fr::from(3u32);
        let a = Fr::from(2u32);
        let b = Fr::from(1u32);
        let c = Fr::from(4u32);
        let d = Fr::from(5u32);
        
        // Compute intermediate values
        let x2 = x * x;                    // 9
        let x3 = x2 * x;                   // 27
        let ax3 = a * x3;                  // 54
        let bx2 = b * x2;                  // 9
        let cx = c * x;                    // 12
        let result = ax3 + bx2 + cx + d;   // 54 + 9 + 12 + 5 = 80
        
        let witness = vec![
            Fr::one(),    // variable 0: constant 1 (public input)
            x,            // variable 1: x = 3 (public input)
            result,       // variable 2: result = 80 (public input)
            a,            // variable 3: a = 2 (private)
            b,            // variable 4: b = 1 (private)
            c,            // variable 5: c = 4 (private)
            d,            // variable 6: d = 5 (private)
            x2,           // variable 7: x2 = 9 (intermediate)
            x3,           // variable 8: x3 = 27 (intermediate)
            ax3,          // variable 9: ax3 = 54 (intermediate)
            bx2,          // variable 10: bx2 = 9 (intermediate)
            cx,           // variable 11: cx = 12 (intermediate)
        ];

        // Public inputs: [1, x, result] - we know x and the final result, but not the coefficients
        let public_inputs = vec![Fr::one(), x, result];

        let proof = prove::<Bn254>(&qap, &witness, &setup);

        assert!(verify::<Bn254>(&qap, &public_inputs, &proof, &setup), "Complex QAP proof verification failed");

        let wrong_public_inputs = vec![Fr::one(), x, result + Fr::one()];
        assert!(!verify::<Bn254>(&qap, &wrong_public_inputs, &proof, &setup), "Proof verification should fail with wrong public input");
    }

    #[test]
    fn test_complex_qap_different_witness() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_complex_qap::<Fr>();
        let setup = trusted_setup::<Bn254>(&qap);

        // Use the same QAP structure but with different polynomial coefficients
        // This time: x=2, a=1, b=3, c=2, d=7
        // Expected: result = 1*8 + 3*4 + 2*2 + 7 = 8 + 12 + 4 + 7 = 31
        
        let x = Fr::from(2u32);
        let a = Fr::from(1u32);
        let b = Fr::from(3u32);
        let c = Fr::from(2u32);
        let d = Fr::from(7u32);
        
        // Compute intermediate values
        let x2 = x * x;                    // 4
        let x3 = x2 * x;                   // 8
        let ax3 = a * x3;                  // 8
        let bx2 = b * x2;                  // 12
        let cx = c * x;                    // 4
        let result = ax3 + bx2 + cx + d;   // 8 + 12 + 4 + 7 = 31
        
        let witness = vec![
            Fr::one(),    // variable 0: constant 1 (public input)
            x,            // variable 1: x = 2 (public input)
            result,       // variable 2: result = 31 (public input)
            a,            // variable 3: a = 1 (private)
            b,            // variable 4: b = 3 (private)
            c,            // variable 5: c = 2 (private)
            d,            // variable 6: d = 7 (private)
            x2,           // variable 7: x2 = 4 (intermediate)
            x3,           // variable 8: x3 = 8 (intermediate)
            ax3,          // variable 9: ax3 = 8 (intermediate)
            bx2,          // variable 10: bx2 = 12 (intermediate)
            cx,           // variable 11: cx = 4 (intermediate)
        ];

        // Public inputs: [1, x, result] - different values than the previous test
        let public_inputs = vec![Fr::one(), x, result];

        let proof = prove::<Bn254>(&qap, &witness, &setup);

        assert!(verify::<Bn254>(&qap, &public_inputs, &proof, &setup), "Complex QAP proof verification failed with different witness");

        // Test with wrong result to ensure security
        let wrong_public_inputs = vec![Fr::one(), x, Fr::from(999u32)];
        assert!(!verify::<Bn254>(&qap, &wrong_public_inputs, &proof, &setup), "Complex QAP should reject wrong result");
    }

    #[test]
    fn test_complex_qap_with_different_curve() {
        use ark_bls12_381::Fr;
        use ark_ff::One;
        
        let qap = create_complex_qap::<Fr>();
        let setup = trusted_setup::<Bls12_381>(&qap);

        // Use the same QAP structure but with different polynomial coefficients
        // This time: x=2, a=1, b=3, c=2, d=7
        // Expected: result = 1*8 + 3*4 + 2*2 + 7 = 8 + 12 + 4 + 7 = 31
        
        let x = Fr::from(2u32);
        let a = Fr::from(1u32);
        let b = Fr::from(3u32);
        let c = Fr::from(2u32);
        let d = Fr::from(7u32);
        
        // Compute intermediate values
        let x2 = x * x;                    // 4
        let x3 = x2 * x;                   // 8
        let ax3 = a * x3;                  // 8
        let bx2 = b * x2;                  // 12
        let cx = c * x;                    // 4
        let result = ax3 + bx2 + cx + d;   // 8 + 12 + 4 + 7 = 31
        
        let witness = vec![
            Fr::one(),    // variable 0: constant 1 (public input)
            x,            // variable 1: x = 2 (public input)
            result,       // variable 2: result = 31 (public input)
            a,            // variable 3: a = 1 (private)
            b,            // variable 4: b = 3 (private)
            c,            // variable 5: c = 2 (private)
            d,            // variable 6: d = 7 (private)
            x2,           // variable 7: x2 = 4 (intermediate)
            x3,           // variable 8: x3 = 8 (intermediate)
            ax3,          // variable 9: ax3 = 8 (intermediate)
            bx2,          // variable 10: bx2 = 12 (intermediate)
            cx,           // variable 11: cx = 4 (intermediate)
        ];

        // Public inputs: [1, x, result] - different values than the previous test
        let public_inputs = vec![Fr::one(), x, result];

        let proof = prove::<Bls12_381>(&qap, &witness, &setup);

        assert!(verify::<Bls12_381>(&qap, &public_inputs, &proof, &setup), "Complex QAP proof verification failed with different witness");

        // Test with wrong result to ensure security
        let wrong_public_inputs = vec![Fr::one(), x, Fr::from(999u32)];
        assert!(!verify::<Bls12_381>(&qap, &wrong_public_inputs, &proof, &setup), "Complex QAP should reject wrong result");
    }

}
