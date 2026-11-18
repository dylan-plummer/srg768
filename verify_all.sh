#!/bin/sh
# Verify every graph file in graphs/ on all pairs (A^2 identity), with the independent numpy verifier.
cd "$(dirname "$0")/code" || exit 1
PY=${PYTHON:-python3}
$PY -B verify_srg.py 768,312,136,120 \
  ../graphs/srg768_312_136_120_threeorbit.g6 \
  ../graphs/srg768_312_136_120_cayley.g6 \
  ../graphs/srg768_312_136_120_norm_threeorbit.g6 \
  ../graphs/srg768_312_136_120_norm_cayley.g6 \
  ../graphs/srg768_312_136_120_family_sample20.g6 &&
$PY -B verify_srg.py 768,455,262,280 ../graphs/srg768_455_262_280_threeorbit_complement.g6
