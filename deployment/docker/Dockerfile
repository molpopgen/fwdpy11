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
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf build

RUN python3 -m pip install --no-cache-dir . \
  && cd / \
  && rm -rf /app
