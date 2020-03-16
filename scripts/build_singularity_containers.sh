#!/usr/bin/env bash

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
BASEDIR=$(dirname $SCRIPTPATH)

CONTAINERS="artic-ncov2019-medaka artic-ncov2019-nanopolish"

for CONTAINER in ${CONTAINERS}; do 
    sudo singularity build ${BASEDIR}/${CONTAINER}.sif ${BASEDIR}/Singularity.${CONTAINER}
done

for CONTAINER in ${CONTAINERS}; do
    echo "Built container ${BASEDIR}/${CONTAINER}"
done

