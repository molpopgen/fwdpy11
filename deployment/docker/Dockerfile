FROM ghcr.io/astral-sh/uv:python3.13-trixie AS build

# Once this image is built:
# If you have a script $HOME/tmp/foo.py,
# you can run it via
# docker run --rm -v $HOME/tmp:/work:z fwdpy11 sh -c ". /app/venv/bin activate && python3 /work/foo.py"
WORKDIR /app

COPY . /app
ENV PATH=/root/.cargo/bin:$PATH

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -qq -y install \
  libgsl-dev \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf build dist \
  # For some reason, installing rustup with apt and then setting the toolchain
  # results in Corrosion being unable to find a valid toolchain.
  # We should try to fix this later.
  && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y \
  && . "$HOME/.cargo/env" \
  # Pin the rustc toolchain to a specific version.
  # Rust 1.64.0 will change the minimum glibc ABI
  # to a version incompatibly with manylinux_2014,
  # so we need to be careful in general.
  && rustup override set 1.66.1 \
  # Pin cbindgen 
  && cargo install --locked cbindgen@0.24.3 \
  && rustc --version \
  && cbindgen --version \
  && uv build . \
  && uv venv /app/venv \
  && . /app/venv/bin/activate \
  && ls -lhrt dist \
  && uv pip install dist/fwdpy11*.whl

FROM python:3.13-slim-trixie
COPY --from=build /app/venv /app/venv
# Make sure the GSL runtime is present
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -qq -y install libgsl28 && apt clean 
