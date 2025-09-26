use std::ops::{AddAssign, Mul, Sub, SubAssign};

use ark_ff::Field;

use crate::pairing_traits::{BasicPairing, BasicPairingGroup};




pub struct LinearPairing<F: Field> {
    _phantom: std::marker::PhantomData<F>,
}

// Wrapper type to avoid conflicts with the generic CurveGroup implementation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TestGroupElement<F: Field>(F);

impl<F: Field> ark_ff::Zero for TestGroupElement<F> {
    fn zero() -> Self {
        TestGroupElement(F::zero())
    }
    
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<F: Field> std::ops::Add<Self> for TestGroupElement<F> {
    type Output = Self;
    
    fn add(self, other: Self) -> Self::Output {
        TestGroupElement(self.0 + other.0)
    }
}

impl<F: Field> AddAssign<Self> for TestGroupElement<F> {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
    }
}

impl<F: Field> SubAssign<Self> for TestGroupElement<F> {
    fn sub_assign(&mut self, other: Self) {
        self.0 -= other.0;
    }
}

impl<F: Field> Sub<Self> for TestGroupElement<F> {
    type Output = Self;
    
    fn sub(self, other: Self) -> Self::Output {
        TestGroupElement(self.0 - other.0)
    }
}

impl<F: Field> Mul<F> for TestGroupElement<F> {
    type Output = Self;
    
    fn mul(self, scalar: F) -> Self::Output {
        TestGroupElement(self.0 * scalar)
    }
}

impl<F: Field> BasicPairingGroup<F> for TestGroupElement<F> {
    fn generator() -> Self {
        TestGroupElement(F::one())
    }
}

impl<F: Field> BasicPairing for LinearPairing<F> {
    type ScalarField = F;
    type G1 = TestGroupElement<F>;
    type G2 = TestGroupElement<F>;
    type TargetGroup = F;
    
    fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetGroup {
        g1_elem.0 * g2_elem.0
    }
}

