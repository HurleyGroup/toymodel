#!/bin/bash

# Params of interest [immutable]
export MU=0.1
export R=0.05
export KN=400000
export KT=200000
export M=40
export CA=23

# execution flags
export CLEAN=true

# Paths of interest
export MODEL_PATH_BASE="/code/model/"
export POST_PATH_BASE="/code/post"
export IMG_PATH_BASE="/code/images"

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


# generate model
python -u /code/src/toymodel.py -mp $MODEL_PATH_BASE
