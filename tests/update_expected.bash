#!/bin/bash

tests="ont_artic
ont_viridian
illumina_artic
illumina_viridian"

for test_name in ${tests}
do
    cp /work/output/${test_name}_test/${test_name}_summary.tsv \
	    /data/pipelines/ncov2019-artic-nf/tests/${test_name}_expected.tsv
done

