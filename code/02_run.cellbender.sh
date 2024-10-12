#! /bin/bash

obj=/object/01_human_heart/snRNA_seq
mkdir ${obj}/CellBender_data/$1
cd ${obj}/CellBender_data/$1
/software/miniconda3/envs/cellbender/bin/cellbender remove-background \
--input ${obj}/CellRanger_data/$1/outs/raw_feature_bc_matrix.h5 \
--output ${obj}/CellBender_data/$1/output.h5 \
--expected-cells $2 \
--total-droplets-included $3 \
