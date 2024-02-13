#!/bin/bash

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

# execution flags
export DEBUG=true

# Paths of interest
export MODEL_PATH_BASE="/code/model/"
export POST_PATH_BASE="/code/post"
export IMG_PATH_BASE="/code/images"

# Paths to save data [hc == hydrostatic constraint | nhc == no hydrostatic constraint]
STRESS_PATH="$POST_PATH_BASE/${PRESSURE}_${LOAD}_${KS_START}_${KS_RATE_RISE}/"
STRESS_IMG="$IMG_PATH_BASE/${PRESSURE}_${LOAD}_${KS_START}_${KS_RATE_RISE}/"

# Explore parameter space (Load, Configuration Angle)
echo "(Load = $LOAD, Configuration Angle = $CA) [PRESSURE = ${PRESSURE}]"

# clear if it exists, and make new folder in its place
if [ "$DEBUG" = false ]; then
  rm -rf $STRESS_PATH
  rm -rf $STRESS_IMG
  mkdir -p $STRESS_PATH
  mkdir -p $STRESS_IMG
fi


# solve for initial conditions
#python -u /code/src/initCond.py -p $STRESS_PATH

# perform integration on DAE
#python -u /code/src/odesolve.py -mp $MODEL_PATH_BASE -p $STRESS_PATH -i $STRESS_IMG --maintain

# post process results
python -u /code/src/postprocess.py -mp $MODEL_PATH_BASE -p $STRESS_PATH -i $STRESS_IMG --maintain

exit 0
