use ark_ff::Field;
use rand_chacha::ChaCha20Rng;
use rand::{SeedableRng, RngCore};

/// Generate toxic waste parameters for Groth16 setup
/// 
/// Returns five field elements: alpha, beta, gamma, delta, x
/// These are the "toxic waste" parameters that must be securely discarded
/// after the trusted setup ceremony.
fn generate_toxic_waste<F: Field>() -> (F, F, F, F, F) {
    // Use ChaCha20: cryptographically secure, fast, side-channel resistant
    let mut rng = ChaCha20Rng::from_entropy();
    
    let alpha = F::rand(&mut rng);
    let beta = F::rand(&mut rng);
    let gamma = F::rand(&mut rng);
    let delta = F::rand(&mut rng);
    let x = F::rand(&mut rng);
    
    (alpha, beta, gamma, delta, x)
}

pub fn setup_groth16<F: Field, const N: usize>() -> (F, F, F, F, F) {
    generate_toxic_waste::<F>()
}

#[cfg(test)]
mod tests {
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
}

