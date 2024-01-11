#!/bin/bash

# Params of interest [immutable]
export MU=0.1
export R=0.05
export KN=400000
export KT=200000
export KS=300000
export M=40

export INT_TIME=0.08

# Paths of interest
export POST_PATH_BASE="/code/post"
export IMG_PATH_BASE="/code/images"

# Explore parameter space (Load, Configuration Angle)
load_start=50
load_end=200
load_step=25

ca_start=10 # degrees
ca_end=50
ca_step=5

load_current=$load_start
while (( $(bc <<< "$load_current < $load_end") )); do
    ca_current=$ca_start
    while (( $(bc <<< "$ca_current < $ca_end") )); do
        echo "(Load = $load_current, Configuration Angle = $ca_current) [MU = $MU]"

        # Paths to save data
        STRESS_RISE_PATH="$POST_PATH_BASE/rise/${MU}/${load_current}_${ca_current}/"
        STRESS_DROP_PATH="$POST_PATH_BASE/drop/${MU}/${load_current}_${ca_current}/"
        STRESS_RISE_IMG="$IMG_PATH_BASE/rise/${MU}/${load_current}_${ca_current}/"
        STRESS_DROP_IMG="$IMG_PATH_BASE/drop/${MU}/${load_current}_${ca_current}/"

        # clear if it exists, and make new folder in its place
        # rm -rf $STRESS_RISE_PATH
        # rm -rf $STRESS_DROP_PATH
        rm -rf $STRESS_RISE_IMG
        rm -rf $STRESS_DROP_IMG

        # mkdir -p $STRESS_RISE_PATH
        # mkdir -p $STRESS_DROP_PATH
        mkdir -p $STRESS_RISE_IMG
        mkdir -p $STRESS_DROP_IMG

        # # run stress rise case
        # python -u /code/src/toymodelRISE.py -p $STRESS_RISE_PATH -l $load_current -a $ca_current # generate functions
        #
        # integrate the ODE and save solution; if it fails, creates a fail file
        # if python -u /code/src/odesolve.py -p $STRESS_RISE_PATH -a $ca_current -d '1E-8'; then
        #     :
        # else
        #     echo "WARNING: Integration failed! Lower integration time?"
        #     echo "$INT_TIME" > $STRESS_RISE_PATH/fail.txt
        # fi

        # postprocess stress rise
        python -u /code/src/postprocessRISE.py -p $STRESS_RISE_PATH -a $ca_current -i $STRESS_RISE_IMG

        # now run stress drop that follows


        # postprocess stress drop

        break
        ca_current=$(bc <<< "$ca_current + $ca_step")
    done
    load_current=$(bc <<< "$load_current + $load_step")
    break
done
