#!/usr/bin/env bash
set -euo pipefail

echo "== Rebuilding grad_coeffs.dat =="

cat tables/grad_coeffs.dat.part_* > tables/grad_coeffs.dat


echo "== Building =="

mkdir -p build
cd build

cmake ..

make -j$(nproc)

echo "== Done =="
