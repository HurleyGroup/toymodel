#!/bin/bash

# Params of interest [immutable]
export MU=0.1
export R=0.05
export KN=400000
export KT=200000
export PRESSURE=350
export KS=1000
export M=40

# Loading conditions
export CA=25
export LOAD=100
export INT_TIME=0.04

# execution flags
export CLEAN=false

# Paths of interest
export POST_PATH_BASE="/code/post"
export IMG_PATH_BASE="/code/images"

# Paths to save data [hc == hydrostatic constraint | nhc == no hydrostatic constraint]
STRESS_PATH="$POST_PATH_BASE/${PRESSURE}/${LOAD}_${CA}/"
STRESS_IMG="$IMG_PATH_BASE/${PRESSURE}/${LOAD}_${CA}/"

# Explore parameter space (Load, Configuration Angle)
echo "(Load = $LOAD, Configuration Angle = $CA) [PRESSURE = ${PRESSURE}]"

# clear if it exists, and make new folder in its place
if [ "$CLEAN" = true ]; then
  rm -rf $POST_PATH_BASE
  rm -rf $IMG_PATH_BASE
  mkdir -p $POST_PATH_BASE
  mkdir -p $IMG_PATH_BASE

  rm -rf $STRESS_PATH
  rm -rf $STRESS_IMG
  mkdir -p $STRESS_PATH
  mkdir -p $STRESS_IMG

fi

# first solve for initial conditions
python -u /code/src/initCond.py -p $STRESS_PATH -a $CA -hp $PRESSURE

# generate model
python -u /code/src/toymodel.py -p $STRESS_PATH -l $LOAD -a $CA -hp $PRESSURE

# run unconstrained case
python -u /code/src/odesolve.py -p $STRESS_PATH -i $STRESS_IMG -a $CA

# run constrained case
python -u /code/src/odesolve.py -p $STRESS_PATH -i $STRESS_IMG -a $CA --maintain

exit 0


#integrate the ODE and save solution; if it fails, creates a fail file
if python -u /code/src/odesolve.py -p $STRESS_RISE_PATH -a $CA  --maintain; then
    :
else
    echo "WARNING: Integration failed! Lower integration time?"
    echo "$INT_TIME_RISE" > $STRESS_RISE_PATH/fail.txt
fi


# run post processing on the stress rise
python -u /code/src/postprocessRISE.py  -p $STRESS_RISE_PATH -i $STRESS_RISE_IMG -a $CA

exit 0
