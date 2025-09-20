use ark_ff::{Field, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_relations::r1cs::ConstraintMatrices;

/// QAP representation 
#[derive(Debug, Clone)]
pub struct QAP<F: Field> {
    pub num_variables: usize,
    pub num_inputs: usize,
    /// u polynomials (one per variable)
    pub u_polynomials: Vec<DensePolynomial<F>>,
    /// V polynomials (one per variable)
    pub v_polynomials: Vec<DensePolynomial<F>>,
    /// W polynomials (one per variable)
    pub w_polynomials: Vec<DensePolynomial<F>>,
    /// Target polynomial T(x) = (x-1)(x-2)...(x-n)
    pub target_polynomial: DensePolynomial<F>,
}

impl<F: Field> QAP<F> {
    /// Generate the division polynomial H(x) = (A(x) * B(x) - C(x)) / T(x)
    /// where A(x), B(x), C(x) are computed from the witness
    pub fn division_polynomial(&self, witness: &[F]) -> DensePolynomial<F> {
        assert_eq!(witness.len(), self.num_variables, "Witness length must match QAP variables");
        
        // Compute U(x) = sum(witness[i] * u_polynomials[i])
        let mut u_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
        for i in 0..self.num_variables {
            let scaled_u = DensePolynomial::from_coefficients_vec(
                self.u_polynomials[i].coeffs().iter().map(|&c| c * witness[i]).collect()
            );
            u_poly = u_poly + scaled_u;
        }
        
        // Compute V(x) = sum(witness[i] * v_polynomials[i])
        let mut v_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
        for i in 0..self.num_variables {
            let scaled_v = DensePolynomial::from_coefficients_vec(
                self.v_polynomials[i].coeffs().iter().map(|&c| c * witness[i]).collect()
            );
            v_poly = v_poly + scaled_v;
        }
        
        // Compute W(x) = sum(witness[i] * w_polynomials[i])
        let mut w_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
        for i in 0..self.num_variables {
            let scaled_w = DensePolynomial::from_coefficients_vec(
                self.w_polynomials[i].coeffs().iter().map(|&c| c * witness[i]).collect()
            );
            w_poly = w_poly + scaled_w;
        }
        
        // Compute U(x) * V(x) - W(x)
        let uv_poly = u_poly.naive_mul(&v_poly);
        let numerator = &uv_poly - &w_poly;
        
        // For a valid witness, numerator should be divisible by T(x) with no remainder
        // Check this by computing quotient * T(x) and seeing if it equals numerator
        let quotient = &numerator / &self.target_polynomial;
        let reconstructed = quotient.naive_mul(&self.target_polynomial);
        
        // Check if the remainder is zero (valid witness)
        if reconstructed != numerator {
            panic!("Invalid witness");
        }
        
        quotient
    }
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

/// Convert arkworks ConstraintMatrices to QAP
/// 
/// # Example
/// ```rust,ignore
/// use ark_bn254::Fr;
/// use ark_relations::r1cs::ConstraintMatrices;
/// use grath::quadratic_arithmetic_programs::constraint_matrices_to_qap;
/// 
/// // Create constraint matrices for x * y = z
/// let matrices = ConstraintMatrices {
///     num_instance_variables: 1,  // constant 1
///     num_witness_variables: 3,   // x, y, z
///     num_constraints: 1,
///     u: vec![vec![(Fr::from(1u64), 1)]], // x coefficient
///     u_num_non_zero: 1,
///     v: vec![vec![(Fr::from(1u64), 2)]], // y coefficient  
///     v_num_non_zero: 1,
///     w: vec![vec![(Fr::from(1u64), 3)]], // z coefficient
///     w_num_non_zero: 1,
/// };
/// 
/// let qap = constraint_matrices_to_qap(&matrices, 1);
/// ```
pub fn constraint_matrices_to_qap<F: Field>(
    matrices: &ConstraintMatrices<F>,
    num_inputs: usize,
) -> QAP<F> {
    let num_constraints = matrices.num_constraints;
    let num_variables = matrices.num_instance_variables + matrices.num_witness_variables;
    
    let mut u_polynomials = Vec::new();
    let mut v_polynomials = Vec::new();
    let mut w_polynomials = Vec::new();
    
    // For each variable, interpolate its column across all constraints
    for var_idx in 0..num_variables {
        // Extract column var_idx from each matrix
        let u_values: Vec<F> = (0..num_constraints)
            .map(|constraint_idx| {
                // Find entry in sparse matrix representation
                matrices.a[constraint_idx]
                    .iter()
                    .find(|(_, col)| *col == var_idx)
                    .map(|(val, _)| *val)
                    .unwrap_or(F::zero())
            })
            .collect();
            
        let v_values: Vec<F> = (0..num_constraints)
            .map(|constraint_idx| {
                matrices.b[constraint_idx]
                    .iter()
                    .find(|(_, col)| *col == var_idx)
                    .map(|(val, _)| *val)
                    .unwrap_or(F::zero())
            })
            .collect();
            
        let w_values: Vec<F> = (0..num_constraints)
            .map(|constraint_idx| {
                matrices.c[constraint_idx]
                    .iter()
                    .find(|(_, col)| *col == var_idx)
                    .map(|(val, _)| *val)
                    .unwrap_or(F::zero())
            })
            .collect();
        
        // Create evaluation points: (1, value[0]), (2, value[1]), ..., (n, value[n-1])
        let points_u: Vec<(F, F)> = (0..num_constraints)
            .map(|i| (F::from((i + 1) as u64), u_values[i]))
            .collect();
            
        let points_v: Vec<(F, F)> = (0..num_constraints)
            .map(|i| (F::from((i + 1) as u64), v_values[i]))
            .collect();
            
        let points_w: Vec<(F, F)> = (0..num_constraints)
            .map(|i| (F::from((i + 1) as u64), w_values[i]))
            .collect();
        
        // Lagrange interpolation
        u_polynomials.push(lagrange_interpolate(&points_u));
        v_polynomials.push(lagrange_interpolate(&points_v));
        w_polynomials.push(lagrange_interpolate(&points_w));
    }
    
    // Create target polynomial T(x) = (x-1)(x-2)...(x-n)
    let target_polynomial = create_target_polynomial(num_constraints);
    
    QAP {
        num_variables,
        num_inputs,
        u_polynomials,
        v_polynomials,
        w_polynomials,
        target_polynomial,
    }
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
    
    #[test]
    fn test_constraint_matrices_to_qap() {
        // Create a simple constraint system: x * y = z
        // Variables: [1, x, y, z] (indices 0, 1, 2, 3)
        // Constraint: x * y = z
        
        let constraint_matrices = ConstraintMatrices {
            num_instance_variables: 1,  // Just the constant 1
            num_witness_variables: 3,   // x, y, z
            num_constraints: 1,
            a: vec![
                vec![(Fr::from(1u64), 1)], // x coefficient in constraint 0
            ],
            a_num_non_zero: 1,
            b: vec![
                vec![(Fr::from(1u64), 2)], // y coefficient in constraint 0
            ],
            b_num_non_zero: 1,
            c: vec![
                vec![(Fr::from(1u64), 3)], // z coefficient in constraint 0
            ],
            c_num_non_zero: 1,
        };
        
        let qap = constraint_matrices_to_qap(&constraint_matrices, 1);
        
        // Verify QAP structure
        assert_eq!(qap.num_variables, 4);
        assert_eq!(qap.num_inputs, 1);
        assert_eq!(qap.u_polynomials.len(), 4);
        assert_eq!(qap.v_polynomials.len(), 4);
        assert_eq!(qap.w_polynomials.len(), 4);
        
        // u-polynomial should be 1 at x=1 for variable 1 (x)
        assert_eq!(qap.u_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // constant term should be 0
        assert_eq!(qap.u_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(1u64));
        assert_eq!(qap.u_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // y term should be 0
        assert_eq!(qap.u_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // z term should be 0
        // v-polynomial should be 1 at x=1 for variable 2 (y)
        assert_eq!(qap.v_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(0u64));
        assert_eq!(qap.v_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(0u64));
        assert_eq!(qap.v_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(1u64));
        assert_eq!(qap.v_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(0u64));
        // w-polynomial should be 1 at x=1 for variable 3 (z)
        assert_eq!(qap.w_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(0u64));
        assert_eq!(qap.w_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(0u64));
        assert_eq!(qap.w_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(0u64));
        assert_eq!(qap.w_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(1u64));
        
        // Target polynomial should evaluate to 0 at constraint evaluation point (x=1)
        assert_eq!(qap.target_polynomial.evaluate(&Fr::from(1u64)), Fr::zero());
    }

    #[test]
    fn test_constraint_matrices_to_qap_addition() {
        // Create a simple constraint system: z + x = w
        // Variables: [1, x, y, z, w] (indices 0, 1, 2, 3, 4)
        // Constraint: (z + x) * 1 = w
        
        let constraint_matrices = ConstraintMatrices {
            num_instance_variables: 1,  // Just the constant 1
            num_witness_variables: 4,   // x, y, z, w
            num_constraints: 1,
            a: vec![
                vec![(Fr::from(1u64), 1), (Fr::from(1u64), 3)], // x and z coefficients in constraint 0
            ],
            a_num_non_zero: 2,
            b: vec![
                vec![(Fr::from(1u64), 0)], // constant term in constraint 0
            ],
            b_num_non_zero: 1,
            c: vec![
                vec![(Fr::from(1u64), 4)], // w coefficient in constraint 0
            ],
            c_num_non_zero: 1,
        };
        
        let qap = constraint_matrices_to_qap(&constraint_matrices, 1);
        
        // Verify QAP structure
        assert_eq!(qap.num_variables, 5);
        assert_eq!(qap.num_inputs, 1);
        assert_eq!(qap.u_polynomials.len(), 5);
        assert_eq!(qap.v_polynomials.len(), 5);
        assert_eq!(qap.w_polynomials.len(), 5);
        
        // u-polynomial should reflect coefficients for x and z
        assert_eq!(qap.u_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // constant term should be 0
        assert_eq!(qap.u_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // x term should be 1
        assert_eq!(qap.u_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // y term should be 0
        assert_eq!(qap.u_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // z term should be 1
        assert_eq!(qap.u_polynomials[4].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // w term should be 0

        // v-polynomial should be 1 at x=1 for constant term
        assert_eq!(qap.v_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // constant term should be 1
        assert_eq!(qap.v_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // x term should be 0
        assert_eq!(qap.v_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // y term should be 0
        assert_eq!(qap.v_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // z term should be 0
        assert_eq!(qap.v_polynomials[4].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // w term should be 0

        // w-polynomial should be 1 at x=1 for w term
        assert_eq!(qap.w_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // constant term should be 0
        assert_eq!(qap.w_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // x term should be 0
        assert_eq!(qap.w_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // y term should be 0
        assert_eq!(qap.w_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // z term should be 0
        assert_eq!(qap.w_polynomials[4].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // w term should be 1 
    }

    #[test]
    fn test_constraint_matrices_to_qap_two_constraints() {
        // Create a simple constraint system:
        // 1) x * y = z
        // 2) z + x = w
        // Variables: [1, x, y, z, w] (indices 0, 1, 2, 3, 4)
        // Constraints:
        // 1) x * y = z
        // 2) (z + x) * 1 = w
        
        let constraint_matrices = ConstraintMatrices {
            num_instance_variables: 1,  // Just the constant 1
            num_witness_variables: 4,   // x, y, z, w
            num_constraints: 2,
            a: vec![
                vec![(Fr::from(1u64), 1)],             // Constraint 0: x coefficient
                vec![(Fr::from(1u64), 1), (Fr::from(1u64), 3)], // Constraint 1: x and z coefficients
            ],
            a_num_non_zero: 3,
            b: vec![
                vec![(Fr::from(1u64), 2)],             // Constraint 0: y coefficient  
                vec![(Fr::from(1u64), 0)],             // Constraint 1: constant term
            ],
            b_num_non_zero: 2,
            c: vec![
                vec![(Fr::from(1u64), 3)],             // Constraint 0: z coefficient
                vec![(Fr::from(1u64), 4)],             // Constraint 1: w coefficient
            ],
            c_num_non_zero: 2,
        };
        
        let qap = constraint_matrices_to_qap(&constraint_matrices, 1);
        
        // Verify QAP structure
        assert_eq!(qap.num_variables, 5);
        assert_eq!(qap.num_inputs, 1);
        assert_eq!(qap.u_polynomials.len(), 5);
        assert_eq!(qap.v_polynomials.len(), 5);
        assert_eq!(qap.w_polynomials.len(), 5);
        
        // u-polynomial should reflect coefficients for x and z across both constraints
        assert_eq!(qap.u_polynomials[1].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // x term should be 1 at constraint 1
        assert_eq!(qap.u_polynomials[1].evaluate(&Fr::from(2u64)), Fr::from(1u64)); // x term should be 1 at constraint 2
        assert_eq!(qap.u_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // z term should be 0 at constraint 1
        assert_eq!(qap.u_polynomials[3].evaluate(&Fr::from(2u64)), Fr::from(1u64)); // z term should be 1 at constraint 2

        // v-polynomial should reflect coefficients for y and constant term across both constraints
        assert_eq!(qap.v_polynomials[2].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // y term should be 1 at constraint 1
        assert_eq!(qap.v_polynomials[2].evaluate(&Fr::from(2u64)), Fr::from(0u64)); // y term should be 0 at constraint 2
        assert_eq!(qap.v_polynomials[0].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // constant term should be 0 at constraint 1
        assert_eq!(qap.v_polynomials[0].evaluate(&Fr::from(2u64)), Fr::from(1u64)); // constant term should be 1 at constraint 2

        // w-polynomial should reflect coefficients for z and w across both constraints
        assert_eq!(qap.w_polynomials[3].evaluate(&Fr::from(1u64)), Fr::from(1u64)); // z term should be 1 at constraint 1
        assert_eq!(qap.w_polynomials[3].evaluate(&Fr::from(2u64)), Fr::from(0u64)); // z term should be 0 at constraint 2
        assert_eq!(qap.w_polynomials[4].evaluate(&Fr::from(1u64)), Fr::from(0u64)); // w term should be 0 at constraint 1
        assert_eq!(qap.w_polynomials[4].evaluate(&Fr::from(2u64)), Fr::from(1u64)); // w term should be 1 at constraint 2
    }

    #[test]
    fn test_division_polynomial() {
        // Create a simple constraint system: x * y = z
        // Variables: [1, x, y, z] (indices 0, 1, 2, 3)
        // Constraint: x * y = z
        
        let constraint_matrices = ConstraintMatrices {
            num_instance_variables: 1,  // Just the constant 1
            num_witness_variables: 3,   // x, y, z
            num_constraints: 1,
            a: vec![
                vec![(Fr::from(1u64), 1)], // x coefficient in constraint 0
            ],
            a_num_non_zero: 1,
            b: vec![
                vec![(Fr::from(1u64), 2)], // y coefficient in constraint 0
            ],
            b_num_non_zero: 1,
            c: vec![
                vec![(Fr::from(1u64), 3)], // z coefficient in constraint 0
            ],
            c_num_non_zero: 1,
        };
        
        let qap = constraint_matrices_to_qap(&constraint_matrices, 1);
        
        // Test with a valid witness: [1, 3, 4, 12] where 3 * 4 = 12
        let witness = vec![Fr::from(1u64), Fr::from(3u64), Fr::from(4u64), Fr::from(12u64)];
        
        // Generate division polynomial
        let h_poly = qap.division_polynomial(&witness);
        
        // For a valid witness, H(x) * T(x) should equal A(x) * B(x) - C(x)
        // Let's verify this at a few points
        for test_point in [Fr::from(5u64), Fr::from(10u64), Fr::from(17u64)] {
            // Compute A(x), B(x), C(x) at test point
            let mut a_val = Fr::zero();
            let mut b_val = Fr::zero();
            let mut c_val = Fr::zero();
            
            for i in 0..qap.num_variables {
                a_val += witness[i] * qap.u_polynomials[i].evaluate(&test_point);
                b_val += witness[i] * qap.v_polynomials[i].evaluate(&test_point);
                c_val += witness[i] * qap.w_polynomials[i].evaluate(&test_point);
            }
            
            let left_side = a_val * b_val - c_val;
            let right_side = h_poly.evaluate(&test_point) * qap.target_polynomial.evaluate(&test_point);
            
            assert_eq!(left_side, right_side, "H(x) * T(x) should equal A(x) * B(x) - C(x) at point {}", test_point);
        }
    }

    #[test]
    #[should_panic(expected = "Invalid witness")]
    fn test_division_polynomial_invalid_witness() {
        // Create a simple constraint system: x * y = z
        // Variables: [1, x, y, z] (indices 0, 1, 2, 3)
        // Constraint: x * y = z
        
        let constraint_matrices = ConstraintMatrices {
            num_instance_variables: 1,  // Just the constant 1
            num_witness_variables: 3,   // x, y, z
            num_constraints: 1,
            a: vec![
                vec![(Fr::from(1u64), 1)], // x coefficient in constraint 0
            ],
            a_num_non_zero: 1,
            b: vec![
                vec![(Fr::from(1u64), 2)], // y coefficient in constraint 0
            ],
            b_num_non_zero: 1,
            c: vec![
                vec![(Fr::from(1u64), 3)], // z coefficient in constraint 0
            ],
            c_num_non_zero: 1,
        };
        
        let qap = constraint_matrices_to_qap(&constraint_matrices, 1);
        
        // Test with an invalid witness: [1, 3, 4, 13] where 3 * 4 != 13
        let witness = vec![Fr::from(1u64), Fr::from(3u64), Fr::from(4u64), Fr::from(13u64)];
        
        let _h_poly = qap.division_polynomial(&witness);
    }
}