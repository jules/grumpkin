<div align="center">
    <h1>
        Grumpkin curve implementation in Rust
    </h1>

**NOTE: While this repository was made mostly to allow for usage of Grumpkin within Halo2, the authors of [halo2curves](https://github.com/privacy-scaling-explorations/halo2curves) have since added it [here](https://github.com/privacy-scaling-explorations/halo2curves/tree/main/src/grumpkin), and using it directly will likely make things easier!**
</div>

This repository implements the [Grumpkin](https://hackmd.io/@aztec-network/ByzgNxBfd#2-Grumpkin---A-curve-on-top-of-BN-254-for-SNARK-efficient-group-operations) curve for use in Rust, by building off of the code provided by [ZCash](https://github.com/zcash) and [Privacy & Scaling Explorations](https://github.com/privacy-scaling-explorations). It should specifically work within [Halo2](https://github.com/zcash/halo2) and proof systems which make use of its components.

## Assembly

This repository extends the use of hand-optimized assembly found in the [halo2curves](https://github.com/privacy-scaling-explorations/halo2curves) repo. To include it, simply use the `asm` feature when importing this crate.
