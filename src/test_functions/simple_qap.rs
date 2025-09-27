//! A simple QAP for testing purposes.
use ark_ff::Field;

use crate::quadratic_arithmetic_programs::QAP;

/// Create a simple QAP that represents the constraint: x * 1 = 4
/// - Constraint: x * 1 = 4
/// 
/// Variables: [1, x] where:
/// - variable 0 is the constant 1 (public input)
/// - variable 1 is x (witness variable)
pub fn qap<F: Field>() -> QAP<F> {
    use ark_relations::r1cs::ConstraintMatrices;
    use crate::quadratic_arithmetic_programs::constraint_matrices_to_qap;
    
    // Variables: [1, x]
    //           [ 0, 1]
    
    let constraint_matrices = ConstraintMatrices {
        num_instance_variables: 1,  // constant 1 (public input)
        num_witness_variables: 1,   // x (witness variable)
        num_constraints: 1,
        a: vec![
            vec![(F::one(), 1)],        // Constraint 0: x * 1 = 4
        ],
        a_num_non_zero: 1,
        b: vec![
            vec![(F::one(), 0)],        // Constraint 0: x * 1 = 4
        ],
        b_num_non_zero: 1,
        c: vec![
            vec![(F::from(4u32), 0)],   // Constraint 0: x * 1 = 4 (output is constant 4)
        ],
        c_num_non_zero: 1,
    };
    
    constraint_matrices_to_qap(&constraint_matrices, 1)
}

/// Create a satisfying witness for the constraint:
/// - Constraint: x * 1 = 4  -> let x=4
pub fn witness<F: Field>() -> Vec<F> {
    vec![F::one(), F::from(4u32)]
}

/// Public inputs are [1] for variable [1]
pub fn public_input<F: Field>() -> Vec<F> {
    vec![F::one()]
}

/// A public input that does not satisfy the QAP for any of the defined witnesses
pub fn wrong_public_input<F: Field>() -> Vec<F> {
    vec![F::from(2u32)] 
}