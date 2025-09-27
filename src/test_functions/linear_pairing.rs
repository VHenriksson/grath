//! # Linear Pairing Implementation for Testing
//!
//! ## Overview
//! 
//! > **⚠️ WARNING: NOT CRYPTOGRAPHICALLY SECURE ⚠️**
//! >
//! > This implementation is **ONLY FOR TESTING** and provides **ZERO CRYPTOGRAPHIC SECURITY**.
//! > The discrete logarithm problem is trivial in this implementation!
//! 
//! This module provides a 'mock' pairing implementation called `LinearPairing`. This is intended only for 
//! testing purposes, having a setting where the discrete logarithm is not just feasible, but trivial.
//! 
//! Instead of two curve groups, we take two copies of a single field, and let the "pairing" be simply
//! multiplication with one element from one copy and one from the other copy.
//! 
//! # 🚨 CRITICAL SECURITY NOTICE 🚨
//! 
//! **DO NOT USE IN PRODUCTION!** This implementation:
//! - Has trivial discrete logarithm
//! - Provides no cryptographic security
//! - Is designed solely for testing and educational purposes
//!
//! ## Usage
//!
//! ```rust,ignore
//! use crate::test_functions::linear_pairing::{LinearPairing, TestGroupElement};
//! use crate::pairing_traits::BasicPairing;
//! use ark_bn254::Fr;
//!
//! type TestPairing = LinearPairing<Fr>;
//! let g1 = TestGroupElement::generator();
//! let mut g2 = TestGroupElement::generator();
//! g2 *= TestGroupElement::generator();
//! let result = TestPairing::pairing(g1, g2);
//! ```
//!
//! ## Example

use std::ops::{AddAssign, Mul, Sub, SubAssign};

use ark_ff::Field;

use crate::pairing_traits::{BasicPairing, BasicPairingGroup};

/// A linear pairing implementation for testing purposes. It uses a single field as both groups,
/// and the pairing operation is simply field multiplication.
///
/// # ⚠️ SECURITY WARNING ⚠️
/// 
/// **THIS PAIRING IS NOT CRYPTOGRAPHICALLY SECURE!** This implementation is
/// designed exclusively for testing and has trivial discrete logarithm properties.
///
/// # Type Parameters
///
/// * `F` - The finite field type used for both group elements and scalars
/// 
/// # Example
/// 
/// ```rust,ignore
/// use ark_bn254::Fr;
/// 
/// type TestPairing = LinearPairing<Fr>;
/// let g1 = TestGroupElement::generator();
/// let g2 = TestGroupElement::generator();
/// let paired = TestPairing::pairing(g1, g2);
/// ```
pub struct LinearPairing<F: Field> {
    /// PhantomData to maintain the generic type parameter without storing data
    _phantom: std::marker::PhantomData<F>,
}

/// A test group element that wraps a field element.
///
/// # Type Parameters
///
/// * `F` - The finite field type being wrapped as a group element
/// 
/// # Example
///
/// ```rust,ignore
/// let g = TestGroupElement::generator();
/// let scalar = Fr::from(42u32);
/// let scaled = g * scalar;
/// let sum = scaled + TestGroupElement::generator();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TestGroupElement<F: Field>(
    /// The underlying field element representing this group element
    pub F
);

impl<F: Field> ark_ff::Zero for TestGroupElement<F> {
    /// Returns the additive identity (zero element) of the group.
    ///
    /// # Returns
    ///
    /// A `TestGroupElement` wrapping the field's zero element.
    fn zero() -> Self {
        TestGroupElement(F::zero())
    }
    
    /// Checks if this group element is the additive identity.
    ///
    /// # Returns
    ///
    /// `true` if this element is zero, `false` otherwise.
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<F: Field> std::ops::Add<Self> for TestGroupElement<F> {
    type Output = Self;
    
    /// Performs group addition of two elements.
    ///
    /// # Arguments
    ///
    /// * `other` - The other group element to add
    ///
    /// # Returns
    ///
    /// A new `TestGroupElement` representing the sum.
    fn add(self, other: Self) -> Self::Output {
        TestGroupElement(self.0 + other.0)
    }
}

impl<F: Field> AddAssign<Self> for TestGroupElement<F> {
    /// Performs in-place group addition.
    ///
    /// # Arguments
    ///
    /// * `other` - The other group element to add
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
    }
}

impl<F: Field> SubAssign<Self> for TestGroupElement<F> {
    /// Performs in-place group subtraction.
    ///
    /// # Arguments
    ///
    /// * `other` - The other group element to subtract
    fn sub_assign(&mut self, other: Self) {
        self.0 -= other.0;
    }
}

impl<F: Field> Sub<Self> for TestGroupElement<F> {
    type Output = Self;
    
    /// Performs group subtraction of two elements.
    ///
    /// # Arguments
    ///
    /// * `other` - The group element to subtract
    ///
    /// # Returns
    ///
    /// A new `TestGroupElement` representing the difference.
    fn sub(self, other: Self) -> Self::Output {
        TestGroupElement(self.0 - other.0)
    }
}

impl<F: Field> Mul<F> for TestGroupElement<F> {
    type Output = Self;
    
    /// Performs scalar multiplication of the group element.
    ///
    /// # Arguments
    ///
    /// * `scalar` - The field element to multiply by
    ///
    /// # Returns
    ///
    /// A new `TestGroupElement` representing the scalar multiple.
    fn mul(self, scalar: F) -> Self::Output {
        TestGroupElement(self.0 * scalar)
    }
}

impl<F: Field> BasicPairingGroup<F> for TestGroupElement<F> {
    /// Returns a generator element for the group.
    ///
    /// # Returns
    ///
    /// A `TestGroupElement` wrapping the field's multiplicative identity.
    fn generator() -> Self {
        TestGroupElement(F::one())
    }
}

impl<F: Field> BasicPairing for LinearPairing<F> {
    /// The scalar field type used for this pairing.
    type ScalarField = F;
    
    /// The first group type (G1) in the pairing.
    type G1 = TestGroupElement<F>;
    
    /// The second group type (G2) in the pairing.
    type G2 = TestGroupElement<F>;
    
    /// The target group type where pairing results reside.
    type TargetGroup = F;
    
    /// Computes the linear pairing between two group elements.
    /// 
    /// # Returns
    ///
    /// A field element from `F` element representing the pairing result.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let g1 = TestGroupElement(Fr::from(3u32));
    /// let g2 = TestGroupElement(Fr::from(4u32));
    /// let result = LinearPairing::pairing(g1, g2); // Fr::from(12u32)
    /// ```
    fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetGroup {
        g1_elem.0 * g2_elem.0
    }
}
