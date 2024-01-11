#!/bin/bash

echo "Building container...."
export DOCKER_BUILDKIT=1
# docker build -t adyota/toymodel ./docker/

echo "Running container...."
docker run -it --rm --name "toy_model_DEBUG" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
#docker run -d --rm --name "toy_model_DETACHED" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"

#
