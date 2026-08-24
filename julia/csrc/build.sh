#!/usr/bin/env bash
# Builds libbinomial_c_api.{dylib,so} from the actual C++ binomial_tree.cpp
# (not a reimplementation), so the Julia tests compare against the real
# reference implementation.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

ROOT=../..
OUT=libbinomial_c_api.dylib
if [[ "$(uname)" != "Darwin" ]]; then
    OUT=libbinomial_c_api.so
fi

c++ -std=c++17 -O2 -shared -fPIC \
    -I"$ROOT/include" \
    binomial_c_api.cpp \
    "$ROOT/src/binomial_tree.cpp" \
    "$ROOT/src/black_scholes.cpp" \
    "$ROOT/src/normal_dist.cpp" \
    -o "$OUT"

echo "built julia/csrc/$OUT"
