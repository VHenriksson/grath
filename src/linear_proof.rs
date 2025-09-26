use ark_ec::pairing::Pairing;
use ark_ff::{Field, UniformRand, Zero, One};
use rand_chacha::ChaCha20Rng;
use rand::{SeedableRng, RngCore};
use ark_poly::Polynomial;
use std::ops::{AddAssign, Mul, Sub, SubAssign, Div};

use crate::quadratic_arithmetic_programs::QAP;
use crate::polynomial_from_exponent_vector::evaluate_polynomial;


pub struct Sigma1<G> {
    pub alpha: G,
    pub beta: G,
    pub delta: G,
    pub x_vec: Vec<G>,
    pub l_vec: Vec<G>,
    pub k_vec: Vec<G>,
    pub t_vec: Vec<G>,
}

pub struct Sigma2<G> {
    pub beta: G,
    pub gamma: G,
    pub delta: G,
    pub x_vec: Vec<G>,
}

/// Groth16 setup parameters
/// 
/// Uses dynamic sizing for flexibility in constraint system sizes
pub struct Groth16SetupParameters<G1, G2> {
    pub sigma1: Sigma1<G1>,
    pub sigma2: Sigma2<G2>,
    // pub alpha: F,
    // pub beta: F,
    // pub gamma: F,
    // pub delta: F,
    // pub x_powers: Vec<F>, // x^0, x^1, ..., x^{N-1}
    // pub l_terms: Vec<F>,  // L terms for input variables
    // pub k_terms: Vec<F>,  // M-L terms for auxiliary variables
    // pub x_powers_times_t_div_by_delta: Vec<F>, // x^i * t(x) / delta for i in 0..N
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

pub fn setup_linear<E: BasicPairing>(qap: &QAP<E::ScalarField>) -> Groth16SetupParameters<E::G1, E::G2> {
    let (alpha, beta, gamma, delta, x) = generate_toxic_waste::<E::ScalarField>();
    setup_linear_with_params::<E>(qap, alpha, beta, gamma, delta, x)
}

pub fn setup_linear_with_params<E: BasicPairing>(
    qap: &QAP<E::ScalarField>,
    alpha: E::ScalarField,
    beta: E::ScalarField,
    gamma: E::ScalarField,
    delta: E::ScalarField,
    x: E::ScalarField
) -> Groth16SetupParameters<E::G1, E::G2> {
    let alpha_1 = E::G1::generator() * alpha;
    let beta_1 = E::G1::generator() * beta;
    let gamma_1 = E::G1::generator() * gamma;
    let delta_1 = E::G1::generator() * delta;
    let alpha_2 = E::G2::generator() * alpha;
    let beta_2 = E::G2::generator() * beta;
    let gamma_2 = E::G2::generator() * gamma;
    let delta_2 = E::G2::generator() * delta;
    
    // Generate powers of x up to the degree needed for the QAP
    let degree = qap.target_polynomial.degree();
    let mut x_powers_1 = Vec::with_capacity(degree + 1);
    let mut x_powers_2 = Vec::with_capacity(degree + 1);
    let mut current_power = E::ScalarField::one();
    for _ in 0..=degree {
        x_powers_1.push(E::G1::generator() * current_power);
        x_powers_2.push(E::G2::generator() * current_power);
        current_power *= x;
    }
    
    // Generate L terms for input variables
    let mut l_terms = vec![E::G1::zero(); qap.num_inputs];
    for i in 0..qap.num_inputs {
        // For Groth16, L terms are typically: (beta * u_i(x) + alpha * v_i(x) + w_i(x)) / gamma
        let u_val = qap.u_polynomials[i].evaluate(&x);
        let v_val = qap.v_polynomials[i].evaluate(&x);
        let w_val = qap.w_polynomials[i].evaluate(&x);
        l_terms[i] = (beta_1 * u_val + alpha_1 * v_val + E::G1::generator() * w_val) / gamma_1;
    }
    
    // Generate K terms for auxiliary variables (witness variables)
    let mut k_terms = vec![E::G1::zero(); qap.num_variables - qap.num_inputs];
    for i in qap.num_inputs..qap.num_variables {
        let u_val = qap.u_polynomials[i].evaluate(&x);
        let v_val = qap.v_polynomials[i].evaluate(&x);
        let w_val = qap.w_polynomials[i].evaluate(&x);
        k_terms[i - qap.num_inputs] = (beta_1 * u_val + alpha_1 * v_val + E::G1::generator() * w_val) / delta_1;
    }
    
    // Generate x^i * t(x) / delta terms
    let mut x_powers_times_t_div_by_delta = Vec::with_capacity(degree);
    if degree >= 1 {
        for i in 0..degree-1 {
            let x_power_i = x_powers_1[i];
            let t_val = qap.target_polynomial.evaluate(&x);
            x_powers_times_t_div_by_delta.push(x_power_i * t_val / delta_1);
        }
    }
    
    Groth16SetupParameters {
        sigma1: Sigma1 {
            alpha: alpha_1,
            beta: beta_1,
            delta: delta_1,
            x_vec: x_powers_1.clone(),
            l_vec: l_terms,
            k_vec: k_terms,
            t_vec: x_powers_times_t_div_by_delta
        },
        sigma2: Sigma2 {
            beta: beta_2,
            gamma: gamma_2,
            delta: delta_2,
            x_vec: x_powers_2
        },
    }
}

pub trait BasicPairingGroup<ScalarField: Field>: 
    Copy + Clone + Zero + 
    AddAssign<Self> + SubAssign<Self> + Sub<Self, Output = Self> + 
    Mul<ScalarField, Output = Self> + Div<Self, Output = Self>
{
    fn generator() -> Self;
    fn inverse(&self) -> Self;
}

pub trait BasicPairing {
    /// The scalar field for the pairing groups
    type ScalarField: Field;
    
    /// First pairing group - must be a module over ScalarField
    type G1: BasicPairingGroup<Self::ScalarField>;
    
    /// Second pairing group - must be a module over ScalarField  
    type G2: BasicPairingGroup<Self::ScalarField>;
    
    /// Target field for pairing results
    type TargetField: Field;
    
    /// Pairing function: takes one element from G1, one from G2, returns element in TargetField
    fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetField;
}

/// NILP Proof structure
#[derive(Debug, Clone)]
struct NILPProof<G1, G2> {
    /// Proof elements for the linear proof
    pub proof_a: G1,  // Evaluation of polynomial A at secret point
    pub proof_b: G2,  // Evaluation of polynomial B at secret point  
    pub proof_c: G1,  // Evaluation of polynomial C at secret point
}

/// NILP Prover - generates a proof that the witness satisfies the QAP
fn prove_linear<E: BasicPairing>(
    qap: &QAP<E::ScalarField>, 
    witness: &[E::ScalarField], 
    setup: &Groth16SetupParameters<E::G1,E::G2>
) -> NILPProof<E::G1,E::G2> {
    let mut rng = &mut ChaCha20Rng::from_entropy();
    let r: E::ScalarField = E::ScalarField::rand(&mut rng);
    let s: E::ScalarField = E::ScalarField::rand(&mut rng);
    prove_linear_with_randomness::<E>(qap, witness, setup, r, s)
}

fn prove_linear_with_randomness<E: BasicPairing>(
    qap: &QAP<E::ScalarField>, 
    witness: &[E::ScalarField], 
    setup: &Groth16SetupParameters<E::G1,E::G2>,
    r: E::ScalarField,
    s: E::ScalarField
) -> NILPProof<E::G1,E::G2> {
    let mut a = setup.sigma1.alpha + setup.sigma1.delta * r;
    let mut b = setup.sigma2.beta + setup.sigma2.delta * s;
    let mut b_g1 = setup.sigma1.beta + setup.sigma1.delta * s;

    for i in 0..witness.len() {
        a += evaluate_polynomial(&qap.u_polynomials[i], &setup.sigma1.x_vec) * witness[i];
        b += evaluate_polynomial(&qap.v_polynomials[i], &setup.sigma2.x_vec) * witness[i];
        b_g1 += evaluate_polynomial(&qap.v_polynomials[i], &setup.sigma1.x_vec) * witness[i];
    }

    let mut c = a * s + b_g1 * r - setup.sigma1.delta * r * s;

    for i in qap.num_inputs..qap.num_variables {
        c += setup.sigma1.k_vec[i - qap.num_inputs] * witness[i];
    }

    c += evaluate_polynomial(&qap.division_polynomial(witness), &setup.sigma1.t_vec);

    NILPProof { proof_a: a, proof_b: b, proof_c: c }
}


fn verify_linear<E: BasicPairing>(
    qap: &QAP<E::ScalarField>,
    public_inputs: &[E::ScalarField],
    proof: &NILPProof<E::G1,E::G2>,
    setup: &Groth16SetupParameters<E::G1,E::G2>
) -> bool {
    assert_eq!(public_inputs.len(), qap.num_inputs, "Public input length must match QAP inputs");
    let mut sum = E::G1::zero();
    for i in 0..qap.num_inputs {
        sum += setup.sigma1.l_vec[i] * public_inputs[i];
    }
    let lhs = E::pairing(proof.proof_a, proof.proof_b);
    let rhs = E::pairing(setup.sigma1.alpha, setup.sigma2.beta) + E::pairing(sum, setup.sigma2.gamma) + E::pairing(proof.proof_c, setup.sigma2.delta);
    lhs == rhs
}


#[cfg(test)]
mod tests {
    use ark_poly::univariate::DensePolynomial;

    use super::*;


    // A fake pairing implementation for testing purposes
    pub struct LinearPairing<F: Field> {
        _phantom: std::marker::PhantomData<F>,
    }

    impl<F: Field> BasicPairingGroup<F> for F {
        fn generator() -> Self {
            F::one()
        }
        
        fn inverse(&self) -> Self {
            self.inverse().unwrap()
        }
    }

    impl<F: Field> BasicPairing for LinearPairing<F> {
        type ScalarField = F;
        type G1 = F;
        type G2 = F;
        type TargetField = F;
        
        fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetField {
            g1_elem * g2_elem
        }
    }


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

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, Fr::zero(), Fr::zero());
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::from(5u32)], &proof, &setup), "Proof verification should fail with wrong public input");
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

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);


        let r = Fr::zero();  // Small non-zero randomness
        let s = Fr::zero();  // Small non-zero randomness

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
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

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        let r = Fr::zero();  // Zero randomness for simplicity
        let s = Fr::zero();  // Zero randomness for simplicity

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");

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

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        // Use non-trivial randomness values
        let r = Fr::from(11u32);  // Non-trivial randomness
        let s = Fr::from(13u32);  // Non-trivial randomness

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
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

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        // Use small but non-trivial randomness values
        let r = Fr::from(13u32);  // Small non-trivial randomness
        let s = Fr::from(17u32);  // Small non-trivial randomness

        // For x=11, our constraint x*1=4 is not satisfied, so we need a different witness
        // Let's use witness [1, 4] but with the QAP evaluated at x=11
        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        let other_proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &other_proof, &setup), "Proof verification failed with all small non-trivial values");
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed with all small non-trivial values");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input");
    }


    #[test]
    fn test_setup_with_true_randomness() {
        use ark_bn254::Fr;
        use ark_ff::{UniformRand, One};
        
        let qap = create_test_qap::<Fr>();
        
        // Use truly random toxic waste parameters as intended in real Groth16
        let (alpha, beta, gamma, delta, x) = generate_toxic_waste::<Fr>();

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        // Use truly random randomness values for the prover
        let mut rng = ChaCha20Rng::from_entropy();
        let r = Fr::rand(&mut rng);
        let s = Fr::rand(&mut rng);

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &[Fr::one(), Fr::from(4u32)], &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::one()], &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &[Fr::from(2u32)], &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
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
    fn test_multi_input_qap_with_fixed_params() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = create_multi_input_qap::<Fr>();

        // Use fixed toxic waste parameters for reproducibility
        let alpha = Fr::from(2u32);   // Small non-trivial value
        let beta = Fr::from(3u32);    // Small non-trivial value
        let gamma = Fr::from(1u32);   // Small non-trivial value
        let delta = Fr::from(2u32);   // Small non-trivial value
        let x = Fr::from(11u32);      // Small non-trivial evaluation point

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

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

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &witness, &setup, Fr::from(13u32), Fr::from(17u32));
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &public_inputs, &proof, &setup), "Multi-input proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &wrong_public_inputs, &proof, &setup), "Multi-input proof verification should fail with wrong public inputs");
    }

    #[test]
    fn test_multi_input_qap() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_multi_input_qap::<Fr>();

        let setup = setup_linear::<LinearPairing<Fr>>(&qap);

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

        let proof = prove_linear::<LinearPairing<Fr>>(&qap, &witness, &setup);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &public_inputs, &proof, &setup), "Multi-input proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &wrong_public_inputs, &proof, &setup), "Multi-input proof verification should fail with wrong public inputs");
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
        

        let setup = setup_linear::<LinearPairing<Fr>>(&qap);

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

        let proof = prove_linear::<LinearPairing<Fr>>(&qap, &witness, &setup);
        let other_proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &witness, &setup, Fr::from(5u32), Fr::from(7u32));
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &public_inputs, &other_proof, &setup), "Proof verification failed");

        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &public_inputs, &proof, &setup), "Ultra-complex QAP proof verification failed");

        let wrong_public_inputs = vec![Fr::one(), x, result + Fr::one()];
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &wrong_public_inputs, &proof, &setup), "Proof verification should fail with wrong public input");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &wrong_public_inputs, &proof, &setup), "Ultra-complex QAP proof verification should fail with wrong public input");
    }

    #[test]
    fn test_complex_qap_different_witness() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = create_complex_qap::<Fr>();
        let setup = setup_linear::<LinearPairing<Fr>>(&qap);

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

        let proof = prove_linear::<LinearPairing<Fr>>(&qap, &witness, &setup);

        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &public_inputs, &proof, &setup), "Complex QAP proof verification failed with different witness");

        // Test with wrong result to ensure security
        let wrong_public_inputs = vec![Fr::one(), x, Fr::from(999u32)];
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &wrong_public_inputs, &proof, &setup), "Complex QAP should reject wrong result");
    }

}