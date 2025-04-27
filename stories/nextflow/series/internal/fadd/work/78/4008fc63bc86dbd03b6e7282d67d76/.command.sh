#!/bin/bash -euo pipefail
gffread \
    annotation.gff \
     \
    --keep-exon-attrs -F -T \
    -o null.gtf

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:PREPARE_GENOME:GFFREAD":
    gffread: $(gffread --version 2>&1)
END_VERSIONS
