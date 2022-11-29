#![cfg_attr(feature = "asm", feature(asm_const))]

#[cfg(feature = "asm")]
mod assembly;

#[macro_use]
mod derive;

mod arithmetic;
mod curve;
mod fq;
mod fr;

#[cfg(test)]
mod tests;

pub use curve::*;
pub use fq::*;
pub use fr::*;

#[cfg(all(feature = "prefetch", target_arch = "x86_64"))]
#[inline(always)]
pub fn prefetch<T>(data: &[T], offset: usize) {
    use core::arch::x86_64::_mm_prefetch;
    unsafe {
        _mm_prefetch(
            data.as_ptr().offset(offset as isize) as *const i8,
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
}
