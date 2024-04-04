#!/bin/bash

# SETUP
set -e -x
rm -rf venv

# INSTALL SYSTEM DEPENDENCIES

yum -y install curl gsl-devel

# INSTALL RUST

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y 
source "$HOME/.cargo/env" 
# Pin the rustc toolchain to a specific version.
# Rust 1.64.0 will change the minimum glibc ABI
# to a version incompatible with manylinux_2014,
# so we need to be careful in general.
rustup override set 1.62.1 
# Pin cbindgen 
cargo install --locked cbindgen@0.24.3 

# Taken from msprime/#2043
# We're running as root in the docker container so git commands issued by
# setuptools_scm will fail without this:
git config --global --add safe.directory /project
