#!/bin/bash

# Params of interest [immutable]
export MU=0.1
export R=0.05
export KN=400000
export KT=200000
export PRESSURE=4550
export M=40

# ks function (should only start changing after ramp finishes)
export KS_START=2000
export KS_RATE_RISE=-400
export KS_TIME_RISE=0.05
export KS_TIME_DROP=0.04
export KS_RATE_DROP=16

# load ramp
export LOAD=0
export LOAD_TIME=.01 # ramp time

# Loading conditions
export CA=23

# execution flags
export CLEAN=false

# Paths of interest
export MODEL_PATH_BASE="/code/model/"
export POST_PATH_BASE="/code/post"
export IMG_PATH_BASE="/code/images"

# Paths to save data [hc == hydrostatic constraint | nhc == no hydrostatic constraint]
STRESS_PATH="$POST_PATH_BASE/${PRESSURE}/${LOAD}_${CA}/"
STRESS_IMG="$IMG_PATH_BASE/${PRESSURE}/${LOAD}_${CA}/"

# Explore parameter space (Load, Configuration Angle)
echo "(Load = $LOAD, Configuration Angle = $CA) [PRESSURE = ${PRESSURE}]"

# clear if it exists, and make new folder in its place
if [ "$CLEAN" = true ]; then
  echo "Cleaning Files."
  rm -rf $POST_PATH_BASE
  rm -rf $IMG_PATH_BASE
  rm -rf $MODEL_PATH_BASE
  mkdir -p $POST_PATH_BASE
  mkdir -p $IMG_PATH_BASE
  mkdir -p $MODEL_PATH_BASE
fi

rm -rf $STRESS_PATH
rm -rf $STRESS_IMG
mkdir -p $STRESS_PATH
mkdir -p $STRESS_IMG


# generate model
# python -u /code/src/toymodel.py -mp $MODEL_PATH_BASE

# solve for initial conditions
python -u /code/src/initCond.py -p $STRESS_PATH

# run unconstrained case
python -u /code/src/odesolve.py -mp $MODEL_PATH_BASE -p $STRESS_PATH -i $STRESS_IMG

# run constrained case
python -u /code/src/odesolve.py -mp $MODEL_PATH_BASE -p $STRESS_PATH -i $STRESS_IMG --maintain


# post process unconstrained case
python -u /code/src/postprocess.py -mp $MODEL_PATH_BASE -p $STRESS_PATH -i $STRESS_IMG

# post process constrained case
python -u /code/src/postprocess.py -mp $MODEL_PATH_BASE -p $STRESS_PATH -i $STRESS_IMG --maintain

exit 0
