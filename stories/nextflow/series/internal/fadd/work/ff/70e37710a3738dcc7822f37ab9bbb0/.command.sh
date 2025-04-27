#!/bin/bash -euo pipefail
[ ! -f  KO_Input_Mock_A_1.fastq.gz ] && ln -sf KO_Input_Mock_A_1.merged.fastq.gz KO_Input_Mock_A_1.fastq.gz
[ ! -f  KO_Input_Mock_A_2.fastq.gz ] && ln -sf KO_Input_Mock_A_2.merged.fastq.gz KO_Input_Mock_A_2.fastq.gz
fastp \
    --in1 KO_Input_Mock_A_1.fastq.gz \
    --in2 KO_Input_Mock_A_2.fastq.gz \
    --out1 KO_Input_Mock_A_1.fastp.fastq.gz \
    --out2 KO_Input_Mock_A_2.fastp.fastq.gz \
    --json KO_Input_Mock_A.fastp.json \
    --html KO_Input_Mock_A.fastp.html \
     \
     \
     \
    --thread 6 \
    --detect_adapter_for_pe \
    --detect_adapter_for_pe --trim_poly_x --trim_poly_g --poly_g_min_len=5 --poly_x_min_len=5 \
    2> >(tee KO_Input_Mock_A.fastp.log >&2)

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTP":
    fastp: $(fastp --version 2>&1 | sed -e "s/fastp //g")
END_VERSIONS
