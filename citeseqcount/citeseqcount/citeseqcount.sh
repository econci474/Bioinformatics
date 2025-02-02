#!/bin/bash

# Exit immediately if there is an error
set -e

SOURCE=/suffolk/RawGenomicsF/ec474/raw/SLX-22953_GEX
HOME=/suffolk/Homes/ec474/scripts
WHITELIST=/suffolk/WorkGenomicsD/ec474/mitox_count_masked/outs/filtered_feature_bc_matrix
OUTPUT=/suffolk/RawGenomicsF/ec474/cite-seq_count5/
echo ENV_SET!

CITE-seq-Count -R1 ${SOURCE}/sample1_R1.fastq.gz \
    -R2 ${SOURCE}/sample1_R2.fastq.gz \
    -t ${HOME}/tags.csv \
    -cbf 1 -cbl 16 -umif 17 -umil 28 \
    -cells 10093 \
    -wl ${WHITELIST}/barcodes.tsv \
    -o ${OUTPUT}
echo CODE_RUN!
