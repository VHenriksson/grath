//! # Basic Pairing Trait Abstractions
//!
//! ## Overview
//!
//! This trait is used to avoid having to implement the full pairing interface for the test pairings. It
//! defines the minimal properties needed for our Groth16 implementation.
//!
//! ## Usage
//!
//! ```rust,ignore
//! use crate::basic_pairing::{BasicPairing, BasicPairingGroup};
//!
//! fn verify<P: BasicPairing>(g1: P::G1, g2: P::G2) -> P::TargetGroup {
//!     P::pairing(g1, g2)
//! }
//! ```

use std::ops::{AddAssign, Mul, Sub, SubAssign};

use ark_ff::{Field, Zero};

/// A trait defining the algebraic structure of groups used in bilinear pairings.
///
/// # Type Parameters
///
/// * `ScalarField` - The finite field over which scalar multiplication is defined
///
/// ```rust,ignore
/// let g = TestGroupElement::generator();
/// let scalar = Fr::from(42u32);
/// let result = g * scalar;
/// ```
pub trait BasicPairingGroup<ScalarField: Field>: 
    Copy + Clone + Zero + std::fmt::Debug +
    AddAssign<Self> + SubAssign<Self> + Sub<Self, Output = Self> + 
    Mul<ScalarField, Output = Self>
{
    /// Returns a generator element for the group.
    ///
    /// # Returns
    ///
    /// A group element that serves as the canonical generator
    fn generator() -> Self;
}

/// The fundamental trait defining bilinear pairing operations.
///
/// ```rust,ignore
/// type TestPairing = LinearPairing<Fr>;
/// let g1 = TestPairing::G1::generator();
/// let g2 = TestPairing::G2::generator();
/// let result = TestPairing::pairing(g1, g2);
/// ```
pub trait BasicPairing {
    /// The scalar field type used for this pairing
    type ScalarField: Field;
    
    /// The first group type (G1) in the pairing
    type G1: BasicPairingGroup<Self::ScalarField>;
    
    /// The second group type (G2) in the pairing
    type G2: BasicPairingGroup<Self::ScalarField>;
    
    /// The target group type where pairing results reside
    type TargetGroup: std::ops::Add<Output = Self::TargetGroup> + PartialEq + std::fmt::Debug;
    
    /// Computes the bilinear pairing between two group elements.
    ///
    /// # Arguments
    ///
    /// * `g1_elem` - An element from the first group (G1)
    /// * `g2_elem` - An element from the second group (G2)
    ///
    /// # Returns
    ///
    /// An element in the target group representing the pairing result
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let g1 = G1::generator() * Fr::from(3u32);
    /// let g2 = G2::generator() * Fr::from(4u32);
    /// let paired = MyPairing::pairing(g1, g2);
    /// ```
    fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetGroup;
}
