FROM ubuntu:focal

# Once this image is built:
# If you have a script $HOME/tmp/foo.py,
# you can run it via
# docker run --rm -v $HOME/tmp:/app fwdpy11 python3 /app/foo.py
WORKDIR /app

COPY . /app

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -qq -y install cmake \
  gcc \
  git \
  g++ \
  libgsl-dev \
  python3 \
  python3-pip \
  curl \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf build \
  && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y \
  && . "$HOME/.cargo/env" \
  # Pin the rustc toolchain to a specific version.
  # Rust 1.64.0 will change the minimum glibc ABI
  # to a version incompatibly with manylinux_2014,
  # so we need to be careful in general.
  && rustup override set 1.62.1 \
  # Pin cbindgen 
  && cargo install cbindgen@0.24.3 \
  && rustc --version \
  && cbindgen --version \
  && ls -lhrt

RUN . "$HOME/.cargo/env" \
  && python3 -m pip install . \
  && cd / \
  && rm -rf /app
