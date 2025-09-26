use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

use crate::quadratic_arithmetic_programs::QAP;



/// Create a QAP that represents the constraint: x * 1 = 4
/// Variables: [1, x] where:
/// - variable 0 is the constant 1 (public input)
/// - variable 1 is x (witness variable)
pub fn qap<F: Field>() -> QAP<F> {
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

pub fn witness<F: Field>() -> Vec<F> {
    vec![F::one(), F::from(4u32)]
}

pub fn public_input<F: Field>() -> Vec<F> {
    vec![F::one()]
}

pub fn wrong_public_input<F: Field>() -> Vec<F> {
    vec![F::from(2u32)] 
}