/// This module is a "hack" to re-export the
/// [demes-forward-capi](https://docs.rs/demes-forward-capi) crate.
///
/// The implementation of the hack is based on
/// [this comment](https://github.com/eqrion/cbindgen/issues/127#issuecomment-362269976).
///
/// Test coverage for this module is upstream.
pub mod demes_capi_bridge;
