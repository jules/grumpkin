#[cfg(feature = "asm")]
use super::assembly::field_arithmetic_asm;
#[cfg(not(feature = "asm"))]
use halo2curves::{field_arithmetic, field_specific};

use crate::arithmetic::{adc, mac, sbb};
use core::convert::TryInto;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use halo2curves::bn256::LegendreSymbol;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// This represents an element of $\mathbb{F}_r$ where
///
/// `p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47`
///
/// is the scalar field of the Grumpkin curve.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Fr` values are always in
// Montgomery form; i.e., Fr(a) = aR mod r, with R = 2^256.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fr(pub(crate) [u64; 4]);

/// Constant representing the modulus
/// q = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
pub const MODULUS: Fr = Fr([
    0x3c208c16d87cfd47,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
]);

/// INV = -(q^{-1} mod 2^64) mod 2^64
const INV: u64 = 0x87d20782e4866389;

/// R = 2^256 mod q
const R: Fr = Fr([
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f,
]);

/// R^2 = 2^512 mod q
const R2: Fr = Fr([
    0xf32cfc5b538afa89,
    0xb5e71911d44501fb,
    0x47ab1eff0a417ff6,
    0x06d89f71cab8351f,
]);

/// R^3 = 2^768 mod q
const R3: Fr = Fr([
    0xb1cd6dafda1530df,
    0x62f210e6a7283db6,
    0xef7f0b0c0ada0afb,
    0x20fd6e902d592544,
]);

/// `GENERATOR = 3 mod r` is a generator of the `r - 1` order multiplicative
/// subgroup, or in other words a primitive root of the field.
const GENERATOR: Fr = Fr::from_raw([0x03, 0x00, 0x00, 0x00]);

/// GENERATOR^t where t * 2^s + 1 = r
/// with t odd. In other words, this
/// is a 2^s root of unity.
const ROOT_OF_UNITY: Fr = Fr([
    0x9e10460b6c3e7ea3,
    0xcbc0b548b438e546,
    0xdc2822db40c0ac2e,
    0x183227397098d014,
]);

const S: u32 = 28;

const MODULUS_STR: &str = "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47";

const TWO_INV: Fr = Fr::from_raw([
    0x9e10460b6c3e7ea4,
    0xcbc0b548b438e546,
    0xdc2822db40c0ac2e,
    0x183227397098d014,
]);

/// 1 / ROOT_OF_UNITY mod r
const ROOT_OF_UNITY_INV: Fr = Fr::zero();

/// GENERATOR^{2^s} where t * 2^s + 1 = r
/// with t odd. In other words, this
/// is a t root of unity.
const DELTA: Fr = Fr::from_raw([
    0xa945f4766d02fffa,
    0x925521d4769f9f48,
    0xaaa1b4e1e5ce6975,
    0x0c93291b6e66ca5a,
]);

/// `ZETA^3 = 1 mod r` where `ZETA^2 != 1 mod r`
const ZETA: Fr = Fr::from_raw([
    0x5763473177fffffeu64,
    0xd4f263f1acdb5c4fu64,
    0x59e26bcea0d48bacu64,
    0x0u64,
]);

use halo2curves::{
    field_common, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output, impl_sum_prod,
};
impl_binops_additive!(Fr, Fr);
impl_binops_multiplicative!(Fr, Fr);
field_common!(
    Fr,
    MODULUS,
    INV,
    MODULUS_STR,
    TWO_INV,
    ROOT_OF_UNITY_INV,
    DELTA,
    ZETA,
    R,
    R2,
    R3
);
impl_sum_prod!(Fr);
#[cfg(not(feature = "asm"))]
field_arithmetic!(Fr, MODULUS, INV, sparse);
#[cfg(feature = "asm")]
field_arithmetic_asm!(Fr, MODULUS, INV);

impl Fr {
    pub const fn size() -> usize {
        32
    }

    pub fn legendre(&self) -> LegendreSymbol {
        // s = self^((modulus - 1) // 2)
        // 0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3
        let s = &[
            0x9e10460b6c3e7ea3u64,
            0xcbc0b548b438e546u64,
            0xdc2822db40c0ac2eu64,
            0x183227397098d014u64,
        ];
        let s = self.pow(s);
        if s == Self::zero() {
            LegendreSymbol::Zero
        } else if s == Self::one() {
            LegendreSymbol::QuadraticResidue
        } else {
            LegendreSymbol::QuadraticNonResidue
        }
    }
}

impl ff::Field for Fr {
    fn random(mut rng: impl RngCore) -> Self {
        Self::from_u512([
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
        ])
    }

    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn double(&self) -> Self {
        self.double()
    }

    #[inline(always)]
    fn square(&self) -> Self {
        self.square()
    }

    /// Computes the square root of this element, if it exists.
    fn sqrt(&self) -> CtOption<Self> {
        let tmp = self.pow(&[
            0x4f082305b61f3f52,
            0x65e05aa45a1c72a3,
            0x6e14116da0605617,
            0x0c19139cb84c680a,
        ]);

        CtOption::new(tmp, !self.ct_eq(&Self::zero()))
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff::helpers::sqrt_ratio_generic(num, div)
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    fn invert(&self) -> CtOption<Self> {
        let tmp = self.pow(&[
            0x3c208c16d87cfd45,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ]);

        CtOption::new(tmp, !self.ct_eq(&Self::zero()))
    }
}

impl ff::PrimeField for Fr {
    type Repr = [u8; 32];

    const NUM_BITS: u32 = 254;
    const CAPACITY: u32 = 253;
    const MODULUS: &'static str = MODULUS_STR;
    const MULTIPLICATIVE_GENERATOR: Self = GENERATOR;
    const ROOT_OF_UNITY: Self = ROOT_OF_UNITY;
    const ROOT_OF_UNITY_INV: Self = ROOT_OF_UNITY_INV;
    const TWO_INV: Self = TWO_INV;
    const DELTA: Self = DELTA;
    const S: u32 = S;

    fn from_repr(repr: Self::Repr) -> CtOption<Self> {
        let mut tmp = Fr([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(repr[0..8].try_into().unwrap());
        tmp.0[1] = u64::from_le_bytes(repr[8..16].try_into().unwrap());
        tmp.0[2] = u64::from_le_bytes(repr[16..24].try_into().unwrap());
        tmp.0[3] = u64::from_le_bytes(repr[24..32].try_into().unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], MODULUS.0[0], 0);
        let (_, borrow) = sbb(tmp.0[1], MODULUS.0[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], MODULUS.0[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], MODULUS.0[3], borrow);

        // If the element is smaller than MODULUS then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    fn to_repr(&self) -> Self::Repr {
        // Turn into canonical form by computing
        // (a.R) / R = a
        #[cfg(feature = "asm")]
        let tmp =
            Self::montgomery_reduce(&[self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0]);

        #[cfg(not(feature = "asm"))]
        let tmp =
            Self::montgomery_reduce(&[self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0]);

        let mut res = [0; 32];
        res[0..8].copy_from_slice(&tmp.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp.0[3].to_le_bytes());

        res
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr()[0] & 1)
    }
}

impl FromUniformBytes<64> for Fr {
    /// Converts a 512-bit little endian integer into
    /// an `Fr` by reducing by the modulus.
    fn from_uniform_bytes(bytes: &[u8; 64]) -> Self {
        Self::from_u512([
            u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
            u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
            u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
            u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
            u64::from_le_bytes(bytes[32..40].try_into().unwrap()),
            u64::from_le_bytes(bytes[40..48].try_into().unwrap()),
            u64::from_le_bytes(bytes[48..56].try_into().unwrap()),
            u64::from_le_bytes(bytes[56..64].try_into().unwrap()),
        ])
    }
}

impl WithSmallOrderMulGroup<3> for Fr {
    const ZETA: Self = ZETA;
}

#[cfg(test)]
mod test {
    use super::*;
    use ff::Field;
    use rand_core::OsRng;

    #[test]
    fn test_sqrt() {
        let v = (Fr::TWO_INV).square().sqrt().unwrap();
        assert!(v == Fr::TWO_INV || (-v) == Fr::TWO_INV);

        for _ in 0..10000 {
            let a = Fr::random(OsRng);
            let mut b = a;
            b = b.square();
            assert_eq!(b.legendre(), LegendreSymbol::QuadraticResidue);

            let b = b.sqrt().unwrap();
            let mut negb = b;
            negb = negb.neg();

            assert!(a == b || a == negb);
        }

        let mut c = Fr::one();
        for _ in 0..10000 {
            let mut b = c;
            b = b.square();
            assert_eq!(b.legendre(), LegendreSymbol::QuadraticResidue);

            b = b.sqrt().unwrap();

            if b != c {
                b = b.neg();
            }

            assert_eq!(b, c);

            c += &Fr::one();
        }
    }

    #[test]
    fn test_field() {
        crate::tests::field::random_field_tests::<Fr>("bn256 scalar".to_string());
    }

    #[test]
    fn test_delta() {
        assert_eq!(Fr::DELTA, GENERATOR.pow(&[1u64 << Fr::S, 0, 0, 0]));
        assert_eq!(
            Fr::DELTA,
            Fr::MULTIPLICATIVE_GENERATOR.pow(&[1u64 << Fr::S, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_u512() {
        assert_eq!(
            Fr::from_raw([
                0x1f8905a172affa8a,
                0xde45ad177dcf3306,
                0xaaa7987907d73ae2,
                0x24d349431d468e30,
            ]),
            Fr::from_u512([
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa
            ])
        );
    }
}
