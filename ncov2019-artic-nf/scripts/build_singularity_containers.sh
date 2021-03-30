#!/usr/bin/env bash

set -euo pipefail

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
BASEDIR=$(dirname $SCRIPTPATH)

WORKFLOWS="illumina nanopore"

for WORKFLOW in ${WORKFLOWS}; do 
    sudo singularity build ${BASEDIR}/artic-ncov2019-${WORKFLOW}.sif ${BASEDIR}/environments/${WORKFLOW}/Singularity
done

for WORKFLOW in ${WORKFLOWS}; do
    echo "Built container ${BASEDIR}/artic-ncov2019-${WORKFLOW}.sif"
done

