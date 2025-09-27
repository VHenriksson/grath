//! A test QAP which is more complex than the simpler examples. Hopefully, it is complex enough
//! to represent a more realistic scenario, although still very simple compared to real-world applications.

use ark_ff::Field;

use crate::quadratic_arithmetic_programs::QAP;



/// Create a complex QAP for testing:
/// 
/// - Constraint 1: x * x = x2          (x^2)
/// - Constraint 2: x2 * x = x3         (x^3) 
/// - Constraint 3: a * x3 = ax3        (a*x^3)
/// - Constraint 4: b * x2 = bx2        (b*x^2)
/// - Constraint 5: c * x = cx          (c*x)
/// - Constraint 6: (ax3+bx2+cx+d) * 1 = result (final sum)
///
/// Variables: \[1, x, result, a, b, c, d, x2, x3, ax3, bx2, cx\] (12 variables)
/// Public inputs: constant 1, x, and result (first 3 variables as required)
pub fn qap<F: Field>() -> QAP<F> {
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

/// Create a satisfying witness for polynomial: result = a*x^3 + b*x^2 + c*x + d
///
/// Let's use: x=3, a=2, b=1, c=4, d=5
/// 
/// Expected: result = 2*27 + 1*9 + 4*3 + 5 = 54 + 9 + 12 + 5 = 80
pub fn witness<F: Field>() -> Vec<F> {

    let x = F::from(3u32);
    let a = F::from(2u32);
    let b = F::from(1u32);
    let c = F::from(4u32);
    let d = F::from(5u32);

    // Compute intermediate values
    let x2 = x * x;                    // 9
    let x3 = x2 * x;                   // 27
    let ax3 = a * x3;                  // 54
    let bx2 = b * x2;                  // 9
    let cx = c * x;                    // 12
    let result = ax3 + bx2 + cx + d;   // 54 + 9 + 12 + 5 = 80

    vec![
        F::one(),     // variable 0: constant 1 (public input)
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
    ]
}

/// Public inputs: [1, x, result]
pub fn public_input<F: Field>() -> Vec<F> {
    let x = F::from(3u32);
    let result = F::from(80u32);
    vec![F::one(), x, result]
}

/// Alternative witness with different polynomial coefficients
/// 
/// Let's use: x=2, a=1, b=3, c=2, d=7  
/// 
/// Expected: result = 1*8 + 3*4 + 2*2 + 7 = 8 + 12 + 4 + 7 = 31
pub fn witness_alternative<F: Field>() -> Vec<F> {
    let x = F::from(2u32);
    let a = F::from(1u32);
    let b = F::from(3u32);
    let c = F::from(2u32);
    let d = F::from(7u32);

    // Compute intermediate values
    let x2 = x * x;                    // 4
    let x3 = x2 * x;                   // 8
    let ax3 = a * x3;                  // 8
    let bx2 = b * x2;                  // 12
    let cx = c * x;                    // 4
    let result = ax3 + bx2 + cx + d;   // 8 + 12 + 4 + 7 = 31

    vec![
        F::one(),     // variable 0: constant 1 (public input)
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
    ]
}

/// Public inputs for the alternative witness: [1, x, result]
pub fn public_input_alternative<F: Field>() -> Vec<F> {
    let x = F::from(2u32);
    let result = F::from(31u32);
    vec![F::one(), x, result]
}

/// A public input that does not satisfy the QAP for any of the defined witnesses
pub fn wrong_public_input<F: Field>() -> Vec<F> {
    let x = F::from(3u32);
    let result = F::from(100u32);
    vec![F::one(), x, result + F::one()]
}
