//! A multi-input QAP for testing multiple public inputs.
//!
//! This module provides a QAP with multiple public inputs to test
//! handling of such cases in proof systems.
use ark_ff::Field;

use crate::quadratic_arithmetic_programs::QAP;

/// Create a QAP with multiple public inputs:
///
/// Constraint 1: x * y = z  
/// 
/// Constraint 2: z * 1 = w
///
/// Variables: [1, x, y, z, w] where:
/// - variable 0 is the constant 1 (public input)
/// - variable 1 is x (public input) 
/// - variables 2,3,4 are y,z,w (witness variables)
pub fn qap<F: Field>() -> QAP<F> {
    use ark_relations::r1cs::ConstraintMatrices;
    use crate::quadratic_arithmetic_programs::constraint_matrices_to_qap;
    
    // Variables: [1, x, y, z, w]
    //           [ 0, 1, 2, 3, 4]
    
    let constraint_matrices = ConstraintMatrices {
        num_instance_variables: 2,  // constant 1 and x (public inputs)
        num_witness_variables: 3,   // y, z, w (witness variables)
        num_constraints: 2,
        a: vec![
            vec![(F::one(), 1)],        // Constraint 0: x * y = z
            vec![(F::one(), 3)],        // Constraint 1: z * 1 = w
        ],
        a_num_non_zero: 2,
        b: vec![
            vec![(F::one(), 2)],        // Constraint 0: x * y = z
            vec![(F::one(), 0)],        // Constraint 1: z * 1 = w
        ],
        b_num_non_zero: 2,
        c: vec![
            vec![(F::one(), 3)],        // Constraint 0: x * y = z
            vec![(F::one(), 4)],        // Constraint 1: z * 1 = w
        ],
        c_num_non_zero: 2,
    };
    
    constraint_matrices_to_qap(&constraint_matrices, 2)
}

/// Create a satisfying witness for both constraints:
/// 
/// Constraint 1: x * y = z  -> let x=3, y=4, then z=12
/// 
/// Constraint 2: z * 1 = w  -> z=12, so w=12
/// 
/// Witness: \[1, 3, 4, 12, 12\] for variables \[1, x, y, z, w\]
pub fn witness<F: Field>() -> Vec<F> {
    vec![
        F::one(),         // variable 0: constant 1 (public input)
        F::from(3u32),    // variable 1: x = 3 (public input)
        F::from(4u32),    // variable 2: y = 4 (witness variable)
        F::from(12u32),   // variable 3: z = 12 (x * y = 3 * 4 = 12) (witness variable)
        F::from(12u32),   // variable 4: w = 12 (z * 1 = 12 * 1 = 12) (witness variable)
    ]
}

/// Public inputs: [1, 3] for variables [constant 1, x]
pub fn public_input<F: Field>() -> Vec<F> {
    vec![F::one(), F::from(3u32)]
}

/// Wrong public inputs for testing verification failures
pub fn wrong_public_input<F: Field>() -> Vec<F> {
    vec![F::one(), F::from(4u32)] // Incorrect x value
}