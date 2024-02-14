#!/bin/bash

echo "Building container...."
# export DOCKER_BUILDKIT=1
# docker build -t adyota/toymodel ./docker/

echo "Running container...."
echo "Creating model...."
#docker run -d --rm --name "toy_model_CREATE" -v $PWD:/code adyota/toymodel "/code/scripts/runMODEL.sh"
#docker wait toy_model_CREATE

echo "Running models...."
ALL_LOADS=490
ALL_KS=6000
HYDROSTATIC=(4460 4560 4660)
ALL_KS_RATE=(-1000 1000)

# vary PRESSURE and KS rate
for h in ${HYDROSTATIC[@]};
do
	for ksr in ${ALL_KS_RATE[@]};
	do
    		echo "toy_model_${ksr}_${h}"
		#docker run -e LOAD=${ALL_LOADS[0]} -e PRESSURE=${HYDROSTATIC[1]} -e KS_START=${ALL_KS[0]} -e KS_RATE_RISE=${ksr} -it --rm --name "toy_model_ALL_KS_RATE_${ksr}" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
		#break
		docker run -e LOAD=${ALL_LOADS} -e PRESSURE=${h} -e KS_START=${ALL_KS} -e KS_RATE_RISE=${ksr} \
			-d --rm --name "M_toy_model_${h}_${ksr}" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"

	done
done


#docker run -e LOAD=650 -e PRESSURE=4550 -e KS_START=9000 -e KS_RATE_RISE=0 -it --rm --name "toy_model_DEBUG" -v $PWD:/code adyota/toymodel "/code/scripts/run.sh"
