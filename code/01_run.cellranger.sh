#! /bin/bash

samplelist=(HA-YF1 HA-YF2 HA-YF3 HA-YF4 HA-YM1 HA-YM2 HA-YM3 HA-YM4 HA-OF1 HA-OF2 HA-OF3 HA-OF4 HA-OM1 HA-OM2 HA-OM3 HA-OM4 H-YF1 H-YF2 H-YF3 H-YF4 H-YM1 H-YM2 H-YM3 H-YM4 H-OF1 H-OF2 H-OF3 H-OF5 H-OM1 H-OM2 H-OM3 H-OM4)

for sample in ${samplelist[*]}
do
/software/CellRanger/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger count \
--id=$sample \
--transcriptome=/software/Reference_geneme/hg19-3.0.0.premrna \
--fastqs=/object/02_human_heart/snRNA_seq/Rawdata/$sample \
--sample=$sample \
--localcores=30 \
--chemistry=SC3Pv3
done