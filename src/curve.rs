use crate::Fq;
use crate::Fr;
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, PrimeField};
use group::Curve;
use group::{prime::PrimeCurveAffine, Group as _, GroupEncoding};
use halo2curves::{Coordinates, CurveAffine, CurveAffineExt, CurveExt};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::new_curve_impl;
use halo2curves::{
    batch_add, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    G1Compressed,
    Fq,
    Fr,
    (G1_GENERATOR_X,G1_GENERATOR_Y),
    G1_B,
    "grumpkin_g1",
);

impl G1 {
    fn endomorphism_base(&self) -> Self {
        unimplemented!();
    }
}

impl CurveAffineExt for G1Affine {
    batch_add!();
}

const G1_GENERATOR_X: Fq = Fq::one();
const G1_GENERATOR_Y: Fq = Fq([
    0x11b2dff1448c41d8,
    0x23d3446f21c77dc3,
    0xaa7b8cf435dfafbb,
    0x14b34cf69dc25d68,
]);
const G1_B: Fq = Fq([
    0xdd7056026000005a,
    0x223fa97acb319311,
    0xcc388229877910c0,
    0x34394632b724eaa,
]);

impl group::cofactor::CofactorGroup for G1 {
    type Subgroup = G1;

    fn clear_cofactor(&self) -> Self {
        *self
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        1.into()
    }
}

#[cfg(test)]
mod tests {
    use crate::G1;

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<G1>();
    }
}
