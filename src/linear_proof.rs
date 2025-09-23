use ark_ff::Field;
use rand_chacha::ChaCha20Rng;
use rand::{SeedableRng, RngCore};
use ark_poly::Polynomial;

use crate::quadratic_arithmetic_programs::QAP;

/// Groth16 setup parameters
/// 
/// Uses dynamic sizing for flexibility in constraint system sizes
pub struct Groth16SetupParameters<F: Field> {
    pub alpha: F,
    pub beta: F,
    pub gamma: F,
    pub delta: F,
    pub x_powers: Vec<F>, // x^0, x^1, ..., x^{N-1}
    pub l_terms: Vec<F>,  // L terms for input variables
    pub k_terms: Vec<F>,  // M-L terms for auxiliary variables
    pub x_powers_times_t_div_by_delta: Vec<F>, // x^i * t(x) / delta for i in 0..N
}

/// Generate toxic waste parameters for Groth16 setup
/// 
/// Returns five field elements: alpha, beta, gamma, delta, x
/// These are the "toxic waste" parameters that must be securely discarded
/// after the trusted setup ceremony.
pub fn generate_toxic_waste<F: Field>() -> (F, F, F, F, F) {
    // Use ChaCha20: cryptographically secure, fast, side-channel resistant
    let mut rng = ChaCha20Rng::from_entropy();
    
    let alpha = F::rand(&mut rng);
    let beta = F::rand(&mut rng);
    let gamma = F::rand(&mut rng);
    let delta = F::rand(&mut rng);
    let x = F::rand(&mut rng);
    
    (alpha, beta, gamma, delta, x)
}

pub fn setup_linear<F: Field>(qap: &QAP<F>) -> Groth16SetupParameters<F> {
    let (alpha, beta, gamma, delta, x) = generate_toxic_waste::<F>();
    setup_linear_with_params(qap, alpha, beta, gamma, delta, x)
}

pub fn setup_linear_with_params<F: Field>(
    qap: &QAP<F>, 
    alpha: F, 
    beta: F, 
    gamma: F, 
    delta: F, 
    x: F
) -> Groth16SetupParameters<F> {
    
    // Generate powers of x up to the degree needed for the QAP
    let degree = qap.target_polynomial.degree();
    let mut x_powers = Vec::with_capacity(degree + 1);
    let mut current_power = F::one();
    for _ in 0..=degree {
        x_powers.push(current_power);
        current_power *= x;
    }
    
    // Generate L terms for input variables
    let mut l_terms = vec![F::zero(); qap.num_inputs];
    for i in 0..qap.num_inputs {
        // For Groth16, L terms are typically: (beta * u_i(x) + alpha * v_i(x) + w_i(x)) / gamma
        let u_val = qap.u_polynomials[i].evaluate(&x);
        let v_val = qap.v_polynomials[i].evaluate(&x);
        let w_val = qap.w_polynomials[i].evaluate(&x);
        l_terms[i] = (beta * u_val + alpha * v_val + w_val) / gamma;
    }
    
    // Generate K terms for auxiliary variables (witness variables)
    let mut k_terms = vec![F::zero(); qap.num_variables - qap.num_inputs];
    for i in qap.num_inputs..qap.num_variables {
        let u_val = qap.u_polynomials[i].evaluate(&x);
        let v_val = qap.v_polynomials[i].evaluate(&x);
        let w_val = qap.w_polynomials[i].evaluate(&x);
        k_terms[i - qap.num_inputs] = (beta * u_val + alpha * v_val + w_val) / delta;
    }
    
    // Generate x^i * t(x) / delta terms
    let mut x_powers_times_t_div_by_delta = Vec::with_capacity(degree);
    if degree >= 2 {
        for i in 0..degree-2 {
            let x_power_i = x_powers[i];
            let t_val = qap.target_polynomial.evaluate(&x);
            x_powers_times_t_div_by_delta.push(x_power_i * t_val / delta);
        }
    }
    
    Groth16SetupParameters {
        alpha,
        beta,
        gamma,
        delta,
        x_powers,
        l_terms,
        k_terms,
        x_powers_times_t_div_by_delta,
    }
}

#[cfg(test)]
mod tests {
    use ark_poly::univariate::DensePolynomial;

    use super::*;

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
    fn test_large_field_no_collisions() {
        use ark_bn254::Fr; // Large cryptographic field
        
        let (alpha, beta, gamma, delta, x) = generate_toxic_waste::<Fr>();
        
        // In a cryptographically large field, collisions should be extremely unlikely
        let elements = vec![alpha, beta, gamma, delta, x];
        let unique_count = elements.iter().collect::<std::collections::HashSet<_>>().len();
        
        assert_eq!(unique_count, 5, "Expected no collisions in large field, but got {} unique elements", unique_count);
    }


    /// NILP Proof structure
    #[derive(Debug, Clone)]
    struct NILPProof<F: Field> {
        /// Proof elements for the linear proof
        pub proof_a: F,  // Evaluation of polynomial A at secret point
        pub proof_b: F,  // Evaluation of polynomial B at secret point  
        pub proof_c: F,  // Evaluation of polynomial C at secret point
    }

    /// NILP Prover - generates a proof that the witness satisfies the QAP
    fn prove_linear<F: Field>(
        qap: &QAP<F>, 
        witness: &[F], 
        setup: &Groth16SetupParameters<F>
    ) -> NILPProof<F> {
        let mut rng = &mut ChaCha20Rng::from_entropy();
        let r: F = F::rand(&mut rng);
        let s: F = F::rand(&mut rng);
        prove_linear_with_randomness(qap, witness, setup, r, s)
    }

    /// NILP Prover with explicit randomness parameters for testing
    fn prove_linear_with_randomness<F: Field>(
        qap: &QAP<F>, 
        witness: &[F], 
        setup: &Groth16SetupParameters<F>,
        r: F,
        s: F
    ) -> NILPProof<F> {
        assert_eq!(witness.len(), qap.num_variables, "Witness length must match QAP variables");
        
        let x = setup.x_powers[1]; // The secret x value (x^1)

        // Calculating h(x) = (U(x)*V(x) - W(x)) / t(x) is actually NOT
        // needed, since we only use it in the context of h(x)t(x), which
        // can be computed directly from U(x), V(x), W(x).
        // We keep it like this for now, since we already implemented it, 
        // and it gives us a free check that the witness satisfies the QAP.
        let h: DensePolynomial<F> = qap.division_polynomial(witness);

        // Compute A(x), B(x), C(x) using the witness
        let mut a_val = F::zero();
        let mut b_val = F::zero();
        let mut c_val = F::zero();

        for i in 0..qap.num_variables {
            let witness_i = witness[i];
            let u_i = qap.u_polynomials[i].evaluate(&x);
            let v_i = qap.v_polynomials[i].evaluate(&x);
            a_val += witness_i * u_i;
            b_val += witness_i * v_i;
            if i >= qap.num_inputs {
                // Only include witness variables in C(x)
                let w_i = qap.w_polynomials[i].evaluate(&x);
                c_val += witness_i * (setup.beta * u_i + setup.alpha * v_i + w_i);
            }
        }

        a_val += setup.alpha + r * setup.delta;
        b_val += setup.beta + s * setup.delta;
        c_val /= setup.delta;
        c_val += qap.target_polynomial.evaluate(&x) * h.evaluate(&x) / setup.delta + s * a_val + r * b_val - r * s * setup.delta;

        
        NILPProof {
            proof_a: a_val,
            proof_b: b_val,
            proof_c: c_val,
        }
    }

    /// NILP Verifier - verifies that a proof is valid for given public inputs
    fn verify_linear<F: Field>(
        qap: &QAP<F>,
        public_inputs: &[F],
        proof: &NILPProof<F>,
        setup: &Groth16SetupParameters<F>
    ) -> bool {
        assert_eq!(public_inputs.len(), qap.num_inputs, "Public input length must match QAP inputs");
        let x = setup.x_powers[1];
        let mut sum = F::zero();
        for i in 0..qap.num_inputs {
            let u_i = qap.u_polynomials[i].evaluate(&x);
            let v_i = qap.v_polynomials[i].evaluate(&x);
            let w_i = qap.w_polynomials[i].evaluate(&x);
            sum += public_inputs[i] * (setup.beta * u_i + setup.alpha * v_i + w_i); //In the linear case, we can ignore gamma
        }
        let lhs = proof.proof_a * proof.proof_b;
        let rhs = setup.alpha * setup.beta + sum + setup.delta * proof.proof_c;
        lhs == rhs
    }

    #[test]
    fn test_setup_with_trivial_values() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = create_test_qap::<Fr>();
        
        let alpha = Fr::zero();
        let beta = Fr::zero();
        let gamma = Fr::one(); 
        let delta = Fr::one();
        let x = Fr::from(4u32); // Evaluation point x=4 (satisfies our constraint)

        let setup = setup_linear_with_params(&qap, alpha, beta, gamma, delta, x);

        // Verify x_powers: should be [1, 4] since x = 4
        assert_eq!(setup.x_powers[0], Fr::one());        // x^0 = 1
        assert_eq!(setup.x_powers[1], Fr::from(4u32));   // x^1 = 4
        
        assert_eq!(setup.l_terms.len(), 1);
        
        assert_eq!(setup.k_terms.len(), 1);

        let proof = prove_linear_with_randomness(&qap, &[Fr::one(), Fr::from(4u32)], &setup, Fr::zero(), Fr::zero());
        assert!(verify_linear(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear(&qap, &[Fr::from(5u32)], &proof, &setup), "Proof verification should fail with wrong public input");
    }
    
    #[test]
    fn test_setup_with_nontrivial_alpha_and_beta() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = create_test_qap::<Fr>();
        
        // Use slightly less trivial toxic waste values
        let alpha = Fr::from(2u32);   // Small non-zero value
        let beta = Fr::from(3u32);    // Small non-zero value
        let gamma = Fr::from(1u32);   // Small non-zero value (avoid zero to prevent division issues)
        let delta = Fr::from(1u32);   // Small non-zero value
        let x = Fr::from(4u32);       // Evaluation point x=4

        let setup = setup_linear_with_params(&qap, alpha, beta, gamma, delta, x);


        let r = Fr::zero();  // Small non-zero randomness
        let s = Fr::zero();  // Small non-zero randomness

        let proof = prove_linear_with_randomness(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
    }

    #[test]
    fn test_setup_with_nontrivial_delta_and_gamma() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = create_test_qap::<Fr>();
        
        // Use trivial alpha and beta, but non-trivial gamma and delta
        let alpha = Fr::zero();       // Trivial value
        let beta = Fr::zero();        // Trivial value
        let gamma = Fr::from(5u32);   // Non-trivial value
        let delta = Fr::from(7u32);   // Non-trivial value
        let x = Fr::from(4u32);       // Evaluation point x=4

        let setup = setup_linear_with_params(&qap, alpha, beta, gamma, delta, x);

        let r = Fr::zero();  // Zero randomness for simplicity
        let s = Fr::zero();  // Zero randomness for simplicity

        let proof = prove_linear_with_randomness(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
        
    }

    #[test]
    fn test_setup_with_nontrivial_randomness() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = create_test_qap::<Fr>();
        
        // Use trivial toxic waste parameters
        let alpha = Fr::zero();       // Trivial value
        let beta = Fr::zero();        // Trivial value
        let gamma = Fr::one();        // Trivial value
        let delta = Fr::one();        // Trivial value
        let x = Fr::from(4u32);       // Evaluation point x=4

        let setup = setup_linear_with_params(&qap, alpha, beta, gamma, delta, x);

        // Use non-trivial randomness values
        let r = Fr::from(11u32);  // Non-trivial randomness
        let s = Fr::from(13u32);  // Non-trivial randomness

        let proof = prove_linear_with_randomness(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
    }

    #[test]
    fn test_setup_with_all_small_nontrivial_values() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_test_qap::<Fr>();
        
        // Use small but non-trivial values for all toxic waste parameters
        let alpha = Fr::from(2u32);   // Small non-trivial value
        let beta = Fr::from(3u32);    // Small non-trivial value
        let gamma = Fr::from(7u32);   // Small non-trivial value
        let delta = Fr::from(2u32);   // Small non-trivial value
        let x = Fr::from(11u32);      // Small non-trivial evaluation point

        let setup = setup_linear_with_params(&qap, alpha, beta, gamma, delta, x);

        // Use small but non-trivial randomness values
        let r = Fr::from(13u32);  // Small non-trivial randomness
        let s = Fr::from(17u32);  // Small non-trivial randomness

        // For x=11, our constraint x*1=4 is not satisfied, so we need a different witness
        // Let's use witness [1, 4] but with the QAP evaluated at x=11
        let proof = prove_linear_with_randomness(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed with all small non-trivial values");
        assert!(!verify_linear(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
    }


    #[test]
    fn test_setup_with_true_randomness() {
        use ark_bn254::Fr;
        use ark_ff::{UniformRand, One};
        
        let qap = create_test_qap::<Fr>();
        
        // Use truly random toxic waste parameters as intended in real Groth16
        let (alpha, beta, gamma, delta, x) = generate_toxic_waste::<Fr>();
        

        let setup = setup_linear_with_params(&qap, alpha, beta, gamma, delta, x);

        // Use truly random randomness values for the prover
        let mut rng = ChaCha20Rng::from_entropy();
        let r = Fr::rand(&mut rng);
        let s = Fr::rand(&mut rng);

        let proof = prove_linear_with_randomness(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify_linear(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
    }

    /// Create a QAP with multiple public inputs:
    /// Constraint 1: x * y = z  
    /// Constraint 2: z * 1 = w
    /// Variables: [1, x, y, z, w] where:
    /// - variable 0 is the constant 1 (public input)
    /// - variable 1 is x (public input) 
    /// - variables 2,3,4 are y,z,w (witness variables)
    fn create_multi_input_qap<F: Field>() -> QAP<F> {
        use ark_poly::DenseUVPolynomial;
        
        // Same constraint structure as the single-input version, but now x is also a public input
        // We have 2 constraints at evaluation points x=1 and x=2
        
        // u polynomials: encode left side of multiplication
        let u_poly_0 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // constant 1: never selected for left side
        let u_poly_1 = DensePolynomial::from_coefficients_vec(vec![F::from(2u32), -F::one()]); // x: selected at x=1, not at x=2 -> 2-x
        let u_poly_2 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // y: never selected for left side
        let u_poly_3 = DensePolynomial::from_coefficients_vec(vec![-F::one(), F::one()]); // z: not at x=1, selected at x=2 -> x-1
        let u_poly_4 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // w: never selected for left side

        // v polynomials: encode right side of multiplication
        let v_poly_0 = DensePolynomial::from_coefficients_vec(vec![-F::one(), F::one()]); // constant 1: not at x=1, selected at x=2 -> x-1
        let v_poly_1 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // x: never selected for right side
        let v_poly_2 = DensePolynomial::from_coefficients_vec(vec![F::from(2u32), -F::one()]); // y: selected at x=1, not at x=2 -> 2-x
        let v_poly_3 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // z: never selected for right side
        let v_poly_4 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // w: never selected for right side

        // w polynomials: encode output of constraint
        let w_poly_0 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // constant 1: never is output
        let w_poly_1 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // x: never is output
        let w_poly_2 = DensePolynomial::from_coefficients_vec(vec![F::zero()]); // y: never is output
        let w_poly_3 = DensePolynomial::from_coefficients_vec(vec![F::from(2u32), -F::one()]); // z: output at x=1, not at x=2 -> 2-x
        let w_poly_4 = DensePolynomial::from_coefficients_vec(vec![-F::one(), F::one()]); // w: not at x=1, output at x=2 -> x-1

        // Target polynomial: (x-1)(x-2) - evaluates to 0 when x=1 or x=2
        let target_poly = DensePolynomial::from_coefficients_vec(vec![
            F::from(2u32),  // constant term: 1*2 = 2
            -F::from(3u32), // x coefficient: -(1+2) = -3
            F::one()        // x^2 coefficient: 1
        ]); // (x-1)(x-2) = x^2 - 3x + 2

        QAP {
            num_variables: 5,  // [1, x, y, z, w]
            num_inputs: 2,     // public inputs: the constant 1 and x
            u_polynomials: vec![u_poly_0, u_poly_1, u_poly_2, u_poly_3, u_poly_4],
            v_polynomials: vec![v_poly_0, v_poly_1, v_poly_2, v_poly_3, v_poly_4],
            w_polynomials: vec![w_poly_0, w_poly_1, w_poly_2, w_poly_3, w_poly_4],
            target_polynomial: target_poly,
        }
    }

    #[test]
    fn test_multi_input_qap() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_multi_input_qap::<Fr>();

        let setup = setup_linear(&qap);

        // Create a satisfying witness for both constraints:
        // Constraint 1: x * y = z  -> let x=3, y=4, then z=12
        // Constraint 2: z * 1 = w  -> z=12, so w=12
        // Witness: [1, 3, 4, 12, 12] for variables [1, x, y, z, w]
        let witness = vec![
            Fr::one(),         // variable 0: constant 1 (public input)
            Fr::from(3u32),    // variable 1: x = 3 (public input)
            Fr::from(4u32),    // variable 2: y = 4 (witness variable)
            Fr::from(12u32),   // variable 3: z = 12 (x * y = 3 * 4 = 12) (witness variable)
            Fr::from(12u32),   // variable 4: w = 12 (z * 1 = 12 * 1 = 12) (witness variable)
        ];

        // Public inputs are [1, 3] for variables [1, x]
        let public_inputs = vec![Fr::one(), Fr::from(3u32)];
        let wrong_public_inputs = vec![Fr::one(), Fr::from(4u32)];

        let proof = prove_linear(&qap, &witness, &setup);
        assert!(verify_linear(&qap, &public_inputs, &proof, &setup), "Multi-input proof verification failed");
        assert!(!verify_linear(&qap, &wrong_public_inputs, &proof, &setup), "Multi-input proof verification should fail with wrong public inputs");
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
        

        let setup = setup_linear(&qap);
        
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

        let proof = prove_linear(&qap, &witness, &setup);
        
        assert!(verify_linear(&qap, &public_inputs, &proof, &setup), "Ultra-complex QAP proof verification failed");
        
        let wrong_public_inputs = vec![Fr::one(), x, result + Fr::one()];
        assert!(!verify_linear(&qap, &wrong_public_inputs, &proof, &setup), "Ultra-complex QAP proof verification should fail with wrong public input");
    }

    #[test]
    fn test_complex_qap_different_witness() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_complex_qap::<Fr>();
        let setup = setup_linear(&qap);
        
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

        let proof = prove_linear(&qap, &witness, &setup);
        
        assert!(verify_linear(&qap, &public_inputs, &proof, &setup), "Complex QAP proof verification failed with different witness");
        
        // Test with wrong result to ensure security
        let wrong_public_inputs = vec![Fr::one(), x, Fr::from(999u32)];
        assert!(!verify_linear(&qap, &wrong_public_inputs, &proof, &setup), "Complex QAP should reject wrong result");
    }
}