use std::ops::{AddAssign, Mul, Sub, SubAssign};

use ark_ff::{Field, Zero};



pub trait BasicPairingGroup<ScalarField: Field>: 
    Copy + Clone + Zero + std::fmt::Debug +
    AddAssign<Self> + SubAssign<Self> + Sub<Self, Output = Self> + 
    Mul<ScalarField, Output = Self>
{
    fn generator() -> Self;
}

pub trait BasicPairing {
    /// The scalar field for the pairing groups
    type ScalarField: Field;
    
    /// First pairing group - must be a module over ScalarField
    type G1: BasicPairingGroup<Self::ScalarField>;
    
    /// Second pairing group - must be a module over ScalarField  
    type G2: BasicPairingGroup<Self::ScalarField>;
    
    /// Target field for pairing results
    type TargetGroup: std::ops::Add<Output = Self::TargetGroup> + PartialEq + std::fmt::Debug;
    
    /// Pairing function: takes one element from G1, one from G2, returns element in TargetField
    fn pairing(g1_elem: Self::G1, g2_elem: Self::G2) -> Self::TargetGroup;
}
