use ark_ff::{Field, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};

/// QAP representation 
#[derive(Debug, Clone)]
pub struct QAP<F: Field> {
    pub num_variables: usize,
    pub num_inputs: usize,
    /// A polynomials (one per variable)
    pub a_polynomials: Vec<DensePolynomial<F>>,
    /// B polynomials (one per variable)
    pub b_polynomials: Vec<DensePolynomial<F>>,
    /// C polynomials (one per variable)
    pub c_polynomials: Vec<DensePolynomial<F>>,
    /// Target polynomial T(x) = (x-1)(x-2)...(x-n)
    pub target_polynomial: DensePolynomial<F>,
}
/// Lagrange interpolation - the heart of R1CS to QAP conversion
fn lagrange_interpolate<F: Field>(points: &[(F, F)]) -> DensePolynomial<F> {
    if points.is_empty() {
        return DensePolynomial::from_coefficients_vec(vec![]);
    }
    
    // First, validate and deduplicate points
    let mut unique_points = Vec::new();
    for &(x, y) in points {
        // Check if we already have a point with this x-coordinate
        if let Some(existing_y) = unique_points.iter().find(|(xi, _)| *xi == x).map(|(_, yi)| *yi) {
            if existing_y != y {
                panic!("Invalid interpolation points: duplicate x-coordinate {} with different y-values {} and {}", 
                       x, existing_y, y);
            }
            // If y-values are the same, skip this duplicate
        } else {
            unique_points.push((x, y));
        }
    }
    
    let points = &unique_points;
    let mut result = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    
    for (i, &(xi, yi)) in points.iter().enumerate() {
        // Create Lagrange basis polynomial for point i
        let mut basis = DensePolynomial::from_coefficients_vec(vec![F::one()]);
        
        for (j, &(xj, _)) in points.iter().enumerate() {
            if i != j {
                let denominator = xi - xj;
                // Now we know denominator cannot be zero since we deduplicated
                // basis *= (x - xj) / (xi - xj)
                let linear = DensePolynomial::from_coefficients_vec(vec![-xj, F::one()]);
                basis = basis.naive_mul(&linear);
                // Multiply by 1/(xi - xj)
                let inv_denom = denominator.inverse().unwrap();
                basis = DensePolynomial::from_coefficients_vec(
                    basis.coeffs().iter().map(|&c| c * inv_denom).collect()
                );
            }
        }
        
        // Add yi * basis to result
        let scaled_basis = DensePolynomial::from_coefficients_vec(
            basis.coeffs().iter().map(|&c| c * yi).collect()
        );
        result = result + scaled_basis;
    }
    
    result
}

/// Create target polynomial T(x) = (x-1)(x-2)...(x-n)
fn create_target_polynomial<F: Field>(num_constraints: usize) -> DensePolynomial<F> {
    let mut target = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    
    for i in 1..=num_constraints {
        // Multiply by (x - i)
        let linear = DensePolynomial::from_coefficients_vec(vec![-F::from(i as u64), F::one()]);
        target = target.naive_mul(&linear);
    }
    
    target
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;
    
    #[test]
    fn test_lagrange_interpolation() {
        use rand::Rng;
        
        // Test interpolating through points (1, 1), (2, 4), (3, 9) -> should give x^2
        let points = vec![
            (Fr::from(1u64), Fr::from(1u64)),
            (Fr::from(2u64), Fr::from(4u64)), 
            (Fr::from(3u64), Fr::from(9u64)),
        ];
        
        let poly = lagrange_interpolate(&points);
        
        // Check that it evaluates correctly at the interpolation points
        assert_eq!(poly.evaluate(&Fr::from(1u64)), Fr::from(1u64));
        assert_eq!(poly.evaluate(&Fr::from(2u64)), Fr::from(4u64));
        assert_eq!(poly.evaluate(&Fr::from(3u64)), Fr::from(9u64));
        
        
        // Test at random points - should match x^2
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let random_val = rng.gen_range(1..1000u64);
            let x = Fr::from(random_val);
            let expected = x * x;  // x^2
            let actual = poly.evaluate(&x);
            assert_eq!(actual, expected, "Failed for x = {}", random_val);
        }
        
        // Test with negative field elements
        for _ in 0..5 {
            let random_val = rng.gen_range(1..100u64);
            let neg_x = -Fr::from(random_val);
            let expected = neg_x * neg_x;  // (-x)^2 = x^2
            let actual = poly.evaluate(&neg_x);
            assert_eq!(actual, expected, "Failed for negative x");
        }
    }
    
    #[test]
    fn test_target_polynomial() {
        use rand::Rng;
        
        // T(x) for 3 constraints should be (x-1)(x-2)(x-3)
        let target = create_target_polynomial::<Fr>(3);
        
        // Should evaluate to 0 at 1, 2, 3
        assert_eq!(target.evaluate(&Fr::from(1u64)), Fr::zero());
        assert_eq!(target.evaluate(&Fr::from(2u64)), Fr::zero());
        assert_eq!(target.evaluate(&Fr::from(3u64)), Fr::zero());
        
        // Should be non-zero at other points
        assert_ne!(target.evaluate(&Fr::from(4u64)), Fr::zero());
        assert_ne!(target.evaluate(&Fr::from(0u64)), Fr::zero());
        
        // Test random points - should be non-zero except at roots 1, 2, 3
        let mut rng = rand::thread_rng();
        for _ in 0..20 {
            let random_val = rng.gen_range(4..1000u64);  // Avoid the roots 1, 2, 3
            let x = Fr::from(random_val);
            let result = target.evaluate(&x);
            assert_ne!(result, Fr::zero(), "Target polynomial should not be zero at x = {}", random_val);
        }
        
        // Test with negative values - should also be non-zero
        for _ in 0..10 {
            let random_val = rng.gen_range(1..100u64);
            let neg_x = -Fr::from(random_val);
            let result = target.evaluate(&neg_x);
            assert_ne!(result, Fr::zero(), "Target polynomial should not be zero at negative x");
        }
    }
    
    #[test]
    #[should_panic(expected = "Invalid interpolation points: duplicate x-coordinate")]
    fn test_lagrange_interpolation_duplicate_x_different_y() {
        // This should panic because we have same x-coordinate with different y-values
        let points = vec![
            (Fr::from(1u64), Fr::from(2u64)),  // (1, 2)
            (Fr::from(1u64), Fr::from(3u64)),  // (1, 3) - same x, different y!
            (Fr::from(2u64), Fr::from(4u64)),
        ];
        
        lagrange_interpolate(&points); // Should panic
    }
    
    #[test]
    fn test_lagrange_interpolation_duplicate_x_same_y() {
        // This should work fine because duplicate points with same y-value are redundant but valid
        let points = vec![
            (Fr::from(1u64), Fr::from(1u64)),  // (1, 1)
            (Fr::from(1u64), Fr::from(1u64)),  // (1, 1) - same point, should be fine
            (Fr::from(2u64), Fr::from(4u64)),  // (2, 4)
            (Fr::from(3u64), Fr::from(9u64)),  // (3, 9)
        ];
        
        let poly = lagrange_interpolate(&points);
        
        // Should still behave like x^2
        assert_eq!(poly.evaluate(&Fr::from(1u64)), Fr::from(1u64));
        assert_eq!(poly.evaluate(&Fr::from(2u64)), Fr::from(4u64));
        assert_eq!(poly.evaluate(&Fr::from(3u64)), Fr::from(9u64));
    }

}