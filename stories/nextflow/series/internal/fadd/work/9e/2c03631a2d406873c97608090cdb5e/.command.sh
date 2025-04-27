#!/bin/bash -euo pipefail
samtools faidx sequence.fa
cut -f 1,2 sequence.fa.fai > sequence.fa.sizes

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:PREPARE_GENOME:CUSTOM_GETCHROMSIZES":
    getchromsizes: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
