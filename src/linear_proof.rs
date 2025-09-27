use ark_ff::{Field, UniformRand, Zero, One};
use rand_chacha::ChaCha20Rng;
use rand::{SeedableRng, RngCore};
use ark_poly::Polynomial;
use std::ops::{AddAssign, Mul, Sub, SubAssign};

use crate::basic_pairing::{BasicPairing, BasicPairingGroup};
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
        let scalar_term = (beta * u_val + alpha * v_val + w_val) / gamma;
        l_terms[i] = E::G1::generator() * scalar_term;
    }
    
    // Generate K terms for auxiliary variables (witness variables)
    let mut k_terms = vec![E::G1::zero(); qap.num_variables - qap.num_inputs];
    for i in qap.num_inputs..qap.num_variables {
        let u_val = qap.u_polynomials[i].evaluate(&x);
        let v_val = qap.v_polynomials[i].evaluate(&x);
        let w_val = qap.w_polynomials[i].evaluate(&x);
        let scalar_term = (beta * u_val + alpha * v_val + w_val) / delta;
        k_terms[i - qap.num_inputs] = E::G1::generator() * scalar_term;
    }
    
    // Generate x^i * t(x) / delta terms
    let mut x_powers_times_t_div_by_delta = Vec::with_capacity(degree);
    if degree >= 1 {
        for i in 0..degree-1 {
            let x_power_i = x_powers_1[i];
            let t_val = qap.target_polynomial.evaluate(&x);
            let scalar_term = t_val / delta;
            x_powers_times_t_div_by_delta.push(x_power_i * scalar_term);
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
/// NILP Proof structure
#[derive(Debug, Clone)]
pub struct NILPProof<G1, G2> {
    /// Proof elements for the linear proof
    pub proof_a: G1,  // Evaluation of polynomial A at secret point
    pub proof_b: G2,  // Evaluation of polynomial B at secret point  
    pub proof_c: G1,  // Evaluation of polynomial C at secret point
}

/// NILP Prover - generates a proof that the witness satisfies the QAP
pub fn prove_linear<E: BasicPairing>(
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


pub fn verify_linear<E: BasicPairing>(
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

    use crate::test_functions::{simple_qap, multi_input_qap, complex_qap};
    use crate::test_functions::linear_pairing::LinearPairing;
    use super::*;


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
        
        let qap = simple_qap::qap::<Fr>();
        
        let alpha = Fr::zero();
        let beta = Fr::zero();
        let gamma = Fr::one(); 
        let delta = Fr::one();
        let x = Fr::from(4u32); // Evaluation point x=4 (satisfies our constraint)

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &simple_qap::witness(), &setup, Fr::zero(), Fr::zero());
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::public_input(), &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input");
    }
    
    #[test]
    fn test_setup_with_nontrivial_alpha_and_beta() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = simple_qap::qap::<Fr>();
        
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
        
        let qap = simple_qap::qap::<Fr>();
        
        // Use trivial alpha and beta, but non-trivial gamma and delta
        let alpha = Fr::zero();       // Trivial value
        let beta = Fr::zero();        // Trivial value
        let gamma = Fr::from(5u32);   // Non-trivial value
        let delta = Fr::from(7u32);   // Non-trivial value
        let x = Fr::from(4u32);       // Evaluation point x=4

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        let r = Fr::zero();  // Zero randomness for simplicity
        let s = Fr::zero();  // Zero randomness for simplicity

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &simple_qap::witness(), &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::public_input(), &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input");

    }

    #[test]
    fn test_setup_with_nontrivial_randomness() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = simple_qap::qap::<Fr>();
        
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

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &simple_qap::witness(), &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::public_input(), &proof, &setup), "Proof verification failed with all small non-trivial values");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input");
    }

    #[test]
    fn test_setup_with_all_small_nontrivial_values() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = simple_qap::qap::<Fr>();
        
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
        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &simple_qap::witness(), &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::public_input(), &proof, &setup), "Proof verification failed with all small non-trivial values");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input");
    }


    #[test]
    fn test_setup_with_true_randomness() {
        use ark_bn254::Fr;
        use ark_ff::{UniformRand, One};
        
        let qap = simple_qap::qap::<Fr>();
        
        // Use truly random toxic waste parameters as intended in real Groth16
        let (alpha, beta, gamma, delta, x) = generate_toxic_waste::<Fr>();

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);

        // Use truly random randomness values for the prover
        let mut rng = ChaCha20Rng::from_entropy();
        let r = Fr::rand(&mut rng);
        let s = Fr::rand(&mut rng);

        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &simple_qap::witness(), &setup, r, s);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::public_input(), &proof, &setup), "Proof verification failed with true randomness");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &simple_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input even with true randomness");
    }

    #[test]
    fn test_multi_input_qap_with_fixed_params() {
        use ark_bn254::Fr;
        use ark_ff::{Zero, One};
        
        let qap = multi_input_qap::qap::<Fr>();

        // Use fixed toxic waste parameters for reproducibility
        let alpha = Fr::from(2u32);   // Small non-trivial value
        let beta = Fr::from(3u32);    // Small non-trivial value
        let gamma = Fr::from(1u32);   // Small non-trivial value
        let delta = Fr::from(2u32);   // Small non-trivial value
        let x = Fr::from(11u32);      // Small non-trivial evaluation point

        let setup = setup_linear_with_params::<LinearPairing<Fr>>(&qap, alpha, beta, gamma, delta, x);


        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &multi_input_qap::witness(), &setup, Fr::from(13u32), Fr::from(17u32));
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &multi_input_qap::public_input(), &proof, &setup), "Multi-input proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &multi_input_qap::wrong_public_input(), &proof, &setup), "Multi-input proof verification should fail with wrong public inputs");
    }

    #[test]
    fn test_multi_input_qap() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = multi_input_qap::qap::<Fr>();

        let setup = setup_linear::<LinearPairing<Fr>>(&qap);


        let proof = prove_linear_with_randomness::<LinearPairing<Fr>>(&qap, &multi_input_qap::witness(), &setup, Fr::from(13u32), Fr::from(17u32));
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &multi_input_qap::public_input(), &proof, &setup), "Multi-input proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &multi_input_qap::wrong_public_input(), &proof, &setup), "Multi-input proof verification should fail with wrong public inputs");
    }

    #[test]
    fn test_complex_qap() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = complex_qap::qap::<Fr>();
        

        let setup = setup_linear::<LinearPairing<Fr>>(&qap);

        let proof = prove_linear::<LinearPairing<Fr>>(&qap, &complex_qap::witness(), &setup);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &complex_qap::public_input(), &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &complex_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input");
    }

    #[test]
    fn test_complex_qap_different_witness() {
        use ark_bn254::Fr;
        use ark_ff::One;
        
        let qap = complex_qap::qap::<Fr>();
        let setup = setup_linear::<LinearPairing<Fr>>(&qap);

        let proof = prove_linear::<LinearPairing<Fr>>(&qap, &complex_qap::witness_alternative(), &setup);
        assert!(verify_linear::<LinearPairing<Fr>>(&qap, &complex_qap::public_input_alternative(), &proof, &setup), "Proof verification failed");
        assert!(!verify_linear::<LinearPairing<Fr>>(&qap, &complex_qap::wrong_public_input(), &proof, &setup), "Proof verification should fail with wrong public input");
    }

}