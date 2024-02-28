#!/bin/bash

# paths
export MODEL_PATH_BASE="/code/model/"


# Params of interest [immutable]
export MU=0.1
export R=0.05
export KN=400000
export KT=200000
export M=40
export CA=22

# ks function (should only start changing after ramp finishes)
export KS_TIME_RISE=3.0
export KS_TIME_DROP=0.00
export KS_RATE_DROP=0.00
export NUM_POINTS=10000

# load ramp
export LOAD_TIME=0.0001 # ramp time


# execute
python -u /code/src/collect.py -mp $MODEL_PATH_BASE
