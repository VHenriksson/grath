use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

use crate::quadratic_arithmetic_programs::QAP;


/// Create a QAP with multiple public inputs:
/// Constraint 1: x * y = z  
/// Constraint 2: z * 1 = w
/// Variables: [1, x, y, z, w] where:
/// - variable 0 is the constant 1 (public input)
/// - variable 1 is x (public input) 
/// - variables 2,3,4 are y,z,w (witness variables)
pub fn qap<F: Field>() -> QAP<F> {
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

// Create a satisfying witness for both constraints:
// Constraint 1: x * y = z  -> let x=3, y=4, then z=12
// Constraint 2: z * 1 = w  -> z=12, so w=12
// Witness: [1, 3, 4, 12, 12] for variables [1, x, y, z, w]
pub fn witness<F: Field>() -> Vec<F> {
    vec![
        F::one(),         // variable 0: constant 1 (public input)
        F::from(3u32),    // variable 1: x = 3 (public input)
        F::from(4u32),    // variable 2: y = 4 (witness variable)
        F::from(12u32),   // variable 3: z = 12 (x * y = 3 * 4 = 12) (witness variable)
        F::from(12u32),   // variable 4: w = 12 (z * 1 = 12 * 1 = 12) (witness variable)
    ]
}

// Public inputs are [1, 3] for variables [1, x]
pub fn public_inputs<F: Field>() -> Vec<F> {
    vec![F::one(), F::from(3u32)]
}

pub fn wrong_public_input<F: Field>() -> Vec<F> {
    vec![F::one(), F::from(4u32)] // Incorrect x value
}