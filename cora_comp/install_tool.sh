#!/bin/bash
# CORA-COMP: install geosets into a venv next to this script. Called as `install_tool.sh v1`.
set -e
cd "$(dirname "$0")/.."  # repo root

SUDO=$([ "$(id -u)" = 0 ] || echo sudo)
if command -v apt-get >/dev/null; then
    $SUDO apt-get update
    $SUDO DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake git python3-venv python3-dev
fi

# ponytail: C++ deps (HiGHS, Eigen, Boost, Qhull) are fetched by CMake; preinstall them to speed up the build
python3 -m venv .cora-venv
.cora-venv/bin/pip install --upgrade pip
.cora-venv/bin/pip install .
.cora-venv/bin/python -c "import geosets_py"
