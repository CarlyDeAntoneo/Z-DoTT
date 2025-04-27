#!/bin/bash -euo pipefail
filter_gtf.py \
    --gtf null.gtf \
    --fasta sequence.fa \
    --prefix sequence

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:PREPARE_GENOME:GTF_FILTER":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
