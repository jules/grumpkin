use crate::arithmetic::mul_512;
use crate::Fq;
use crate::Fr;
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, PrimeField};
use group::Curve;
use group::{prime::PrimeCurveAffine, Group as _, GroupEncoding};
use halo2curves::{Coordinates, CurveAffine, CurveAffineExt, CurveExt, FieldExt, Group};
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

const ENDO_G1: [u64; 4] = [
    0x7a7bd9d4391eb18du64,
    0x4ccef014a773d2cfu64,
    0x0000000000000002u64,
    0u64,
];
const ENDO_G2: [u64; 4] = [0xd91d232ec7e0b3d7u64, 0x0000000000000002u64, 0u64, 0u64];
const ENDO_MINUS_B1: [u64; 4] = [0x8211bbeb7d4f1128u64, 0x6f4d8248eeb859fcu64, 0u64, 0u64];
const ENDO_B2: [u64; 4] = [0x89d3256894d213e3u64, 0u64, 0u64, 0u64];
const ENDO_BETA: Fr = Fr::from_raw([
    0xe4bd44e5607cfd48,
    0xc28f069fbb966e3d,
    0x5e6dd9e7e0acccb0,
    0x30644e72e131a029,
]);

trait CurveEndo: CurveExt {
    fn endomorphism_base(&self) -> Self;
    fn endomorphism_scalars(k: &Self::ScalarExt) -> (u128, u128);
}

impl CurveEndo for G1 {
    fn endomorphism_base(&self) -> Self {
        Self {
            x: self.x * Self::Base::ZETA,
            y: -self.y,
            z: self.z,
        }
    }

    fn endomorphism_scalars(k: &Self::ScalarExt) -> (u128, u128) {
        #[cfg(feature = "asm")]
        let input = Fr::montgomery_reduce(&[k.0[0], k.0[1], k.0[2], k.0[3], 0, 0, 0, 0]).0;

        #[cfg(not(feature = "asm"))]
        let input = Fr::montgomery_reduce(k.0[0], k.0[1], k.0[2], k.0[3], 0, 0, 0, 0).0;

        let c1_512 = mul_512(ENDO_G2, input);
        let c2_512 = mul_512(ENDO_G1, input);

        let c1_hi = [c1_512[4], c1_512[5], c1_512[6], c1_512[7]];
        let c2_hi = [c2_512[4], c2_512[5], c2_512[6], c2_512[7]];

        let q1_512 = mul_512(c1_hi, ENDO_MINUS_B1);
        let q2_512 = mul_512(c2_hi, ENDO_B2);

        let q1_lo = Self::ScalarExt::from_raw([q1_512[0], q1_512[1], q1_512[2], q1_512[3]]);
        let q2_lo = Self::ScalarExt::from_raw([q2_512[0], q2_512[1], q2_512[2], q2_512[3]]);

        let k1 = q2_lo - q1_lo;
        let k2 = (k1 * ENDO_BETA) + k;

        (k2.get_lower_128(), k1.get_lower_128())
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        curve::{CurveEndo, ENDO_BETA},
        Fr, G1Affine, G1,
    };
    use ff::Field;
    use halo2curves::CurveExt;
    use rand_core::OsRng;

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<G1>();
    }

    #[test]
    fn test_endo_consistency() {
        let g = G1::generator();
        assert_eq!(g * (-ENDO_BETA), g.endo());
    }

    #[test]
    fn test_endomorphism() {
        use halo2curves::FieldExt;

        let scalar = Fr::random(OsRng);
        let point = G1Affine::random(OsRng);

        let expected = point * scalar;

        let (part1, part2) = G1::endomorphism_scalars(&scalar);

        let k1 = Fr::from_u128(part1);
        let k2 = Fr::from_u128(part2);

        let t1 = point * k1;
        let base = G1::endomorphism_base(&point.into());

        let t2 = base * k2;
        let result = t1 + t2;

        let res_affine: G1Affine = result.into();
        let exp_affine: G1Affine = expected.into();

        assert_eq!(res_affine, exp_affine);
    }
}
