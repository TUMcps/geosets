FROM --platform=linux/amd64 debian:trixie-slim

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    cmake \
    make \
    ninja-build \
    git \
    # Dependencies
    libboost-dev \
    libeigen3-dev \
    libqhull-dev \
    libcdd-dev \
  && rm -rf /var/lib/apt/lists/*

# Build HiGHS from source — Debian's libhighs-dev has broken cmake config
# (DESTDIR leak + relative CMAKE_INSTALL_PREFIX, see bugs.debian.org #1128055)
ARG HIGHS_VERSION=v1.10.0
RUN git clone --depth 1 --branch ${HIGHS_VERSION} https://github.com/ERGO-Code/HiGHS.git /tmp/highs \
 && cmake -S /tmp/highs -B /tmp/highs/build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DHIGHS_BUILD_CLI=OFF \
      -DHIGHS_BUILD_TESTING=OFF \
      -DHIGHS_DEV=OFF \
 && cmake --build /tmp/highs/build --target install \
 && ldconfig \
 && rm -rf /tmp/highs
