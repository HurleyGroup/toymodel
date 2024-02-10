#!/bin/bash

echo "Building container...."
export DOCKER_BUILDKIT=1
# docker build -t adyota/toymodel ./docker/

echo "Running container...."
echo "Creating model...."
docker run -d --rm --name "toy_model_CREATE" -v $PWD:/code adyota/toymodel "/code/scripts/runMODEL.sh"
docker wait toy_model_CREATE

echo "Running models...."
ALL_LOADS=(450 475 490 505 520)
HYDROSTATIC=(1000 4460 4560 4650)
ALL_KS=(3000 6000 9000 12000 15000)
ALL_KS_RATE=(-1000 -500 -100 0 100 500 1000)

# vary LOAD
for l in ${ALL_LOADS[@]};
do
    echo "toy_model_LOAD_${l}"
    docker run -e LOAD=${l} -e PRESSURE=${HYDROSTATIC[2]} -e KS_START=${ALL_KS[2]} -e KS_RATE_RISE=${ALL_KS_RATE[2]} -d --rm --name "toy_model_LOAD_${l}" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
done

# vary PRESSURE
for h in ${HYDROSTATIC[@]};
do
    echo "toy_model_HYDROSTATIC_${h}"
    docker run -e LOAD=${ALL_LOADS[2]} -e PRESSURE=${h} -e KS_START=${ALL_KS[2]} -e KS_RATE_RISE=${ALL_KS_RATE[2]} -d --rm --name "toy_model_HYDROSTATIC_${h}" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
done


# vary KS value
for ks in ${ALL_KS[@]};
do
    echo "toy_model_ALL_KS_${ks}"
    docker run -e LOAD=${ALL_LOADS[2]} -e PRESSURE=${HYDROSTATIC[2]} -e KS_START=${ks} -e KS_RATE_RISE=${ALL_KS_RATE[2]} -d --rm --name "toy_model_ALL_KS_${ks}" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
done


# vary KS rate
for ksr in ${ALL_KS_RATE[@]};
do
    echo "toy_model_ALL_KS_RATE_${ksr}"
    docker run -e LOAD=${ALL_LOADS[2]} -e PRESSURE=${HYDROSTATIC[2]} -e KS_START=${ALL_KS[2]} -e KS_RATE_RISE=${ksr} -d --rm --name "toy_model_ALL_KS_RATE_${ksr}" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
done


# docker run -e LOAD=0 -e PRESSURE=4550 -e KS_START=1000 -e KS_RATE_RISE=50 -it --rm --name "toy_model_DEBUG" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
