//! # Polynomial Evaluation with Exponent Vectors
//!
//! ## Overview
//!
//! A key feature of the Groth16 proof system is the ability to evaluate polynomials without multiplication
//! being available to us. This is achieved by using a precomputed vector of group elements which is delivered
//! from the trusted setup phase. 
//! 
//! Here we implement a function to evaluate a polynomial given such a vector.
//!
//! ## Usage
//!
//! ```rust,ignore
//! use ark_poly::DenseUVPolynomial;
//! use crate::polynomial_from_exponent_vector::evaluate_polynomial;
//!
//! let power_vector = vec![g, x_g, x2_g, x3_g]; // [G, x*G, x²*G, x³*G]
//! let result = evaluate_polynomial(&poly, &power_vector);
//! ```

use ark_ff::{Field, Zero};
use std::ops::{Mul, AddAssign};

/// Evaluate a polynomial using a precomputed power vector.
///
/// # Type Parameters
///
/// * `F` - The finite field type for polynomial coefficients
/// * `G` - The group type for the power vector elements
///
/// # Arguments
///
/// * `poly` - The polynomial to evaluate
/// * `power_vector` - Precomputed powers [1*G, x*G, x²*G, ...]
///
/// # Returns
///
/// The polynomial evaluation result in the group G
///
/// # Example
///
/// ```rust,ignore
/// let poly = DenseUVPolynomial::from_coefficients_vec(vec![Fr::from(7u32), Fr::from(13u32)]);
/// let power_vector = vec![g, x_g]; // [G, x*G]
/// let result = evaluate_polynomial(&poly, &power_vector); // (7 + 13*x)*G
/// ```
pub fn evaluate_polynomial<F: Field, G>(
    poly: &ark_poly::univariate::DensePolynomial<F>, 
    power_vector: &[G]
) -> G 
where
    G: Clone + Zero + AddAssign<G>,
    G: Mul<F, Output = G>,
{
    // Evaluate polynomial p(x) = a₀ + a₁x + a₂x² + ... using powers [1, x, x², ...]
    // Result is a₀*G + a₁*(x*G) + a₂*(x²*G) + ...
    
    let coeffs = &poly.coeffs;  // Access coefficients directly
    let mut result = G::zero();
    
    for (coeff, power_element) in coeffs.iter().zip(power_vector.iter()) {
        result += power_element.clone() * *coeff;
    }
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{G1Projective, Fr};
    use ark_ec::Group;

    #[test]
    fn test_evaluate_polynomial() {
        use ark_poly::DenseUVPolynomial;
        use ark_poly::Polynomial;
        
        // Create a much more complex polynomial p(x) = 7 + 13x - 5x² + 11x³ - 3x⁴ + 2x⁵ + 17x⁶
        let poly = DenseUVPolynomial::from_coefficients_vec(vec![
            Fr::from(7u32),   // constant term
            Fr::from(13u32),  // x term  
            Fr::from(-5i32),  // x² term (negative coefficient)
            Fr::from(11u32),  // x³ term
            Fr::from(-3i32),  // x⁴ term (negative coefficient)
            Fr::from(2u32),   // x⁵ term
            Fr::from(17u32),  // x⁶ term
        ]);


        // Choose a known secret value x = 7
        let x = Fr::from(7u32);
        let g1 = G1Projective::generator();
        
        // Create power vector: [G, x*G, x²*G, ..., x⁶*G] = [G, 7*G, 49*G, 343*G, 2401*G, 16807*G, 117649*G]
        let mut x_power = Fr::from(1u32); // Start with x^0 = 1
        let mut power_vector = Vec::new();
        
        for _ in 0..7 {  // Need 7 powers for degree 6 polynomial
            power_vector.push(g1 * x_power);
            x_power *= x; // Next power: x^1, x^2, x^3, x^4, x^5, x^6
        }
        
        // Method 1: Use our evaluate_polynomial function
        // Should compute: 7*(1*G) + 13*(7*G) - 5*(49*G) + 11*(343*G) - 3*(2401*G) + 2*(16807*G) + 17*(117649*G)
        let result1 = evaluate_polynomial(&poly, &power_vector);
        
        // Method 2: Direct evaluation p(7) * G
        let poly_value = poly.evaluate(&x); 
        // p(7) = 7 + 13*7 - 5*7² + 11*7³ - 3*7⁴ + 2*7⁵ + 17*7⁶
        // Let's verify a few terms manually:
        // 7 + 91 - 5*49 + 11*343 - 3*2401 + 2*16807 + 17*117649
        // = 7 + 91 - 245 + 3773 - 7203 + 33614 + 2000033 = 2030070
        let result2 = g1 * poly_value;
        
        // They should be equal: both represent p(7)*G
        assert_eq!(result1, result2);
        
        // Verify the actual computation manually
        let expected = Fr::from(7u32) + 
                      Fr::from(13u32) * Fr::from(7u32) + 
                      Fr::from(-5i32) * Fr::from(49u32) + 
                      Fr::from(11u32) * Fr::from(343u32) + 
                      Fr::from(-3i32) * Fr::from(2401u32) + 
                      Fr::from(2u32) * Fr::from(16807u32) + 
                      Fr::from(17u32) * Fr::from(117649u32);
        
        assert_eq!(poly_value, expected);
    }
}