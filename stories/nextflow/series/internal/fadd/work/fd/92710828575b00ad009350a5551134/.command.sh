#!/bin/bash -euo pipefail
[ ! -f  WT_Input_Mock_C_1.fastq.gz ] && ln -sf WT_Input_Mock_C_1.merged.fastq.gz WT_Input_Mock_C_1.fastq.gz
[ ! -f  WT_Input_Mock_C_2.fastq.gz ] && ln -sf WT_Input_Mock_C_2.merged.fastq.gz WT_Input_Mock_C_2.fastq.gz
fastp \
    --in1 WT_Input_Mock_C_1.fastq.gz \
    --in2 WT_Input_Mock_C_2.fastq.gz \
    --out1 WT_Input_Mock_C_1.fastp.fastq.gz \
    --out2 WT_Input_Mock_C_2.fastp.fastq.gz \
    --json WT_Input_Mock_C.fastp.json \
    --html WT_Input_Mock_C.fastp.html \
     \
     \
     \
    --thread 6 \
    --detect_adapter_for_pe \
    --detect_adapter_for_pe --trim_poly_x --trim_poly_g --poly_g_min_len=5 --poly_x_min_len=5 \
    2> >(tee WT_Input_Mock_C.fastp.log >&2)

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTP":
    fastp: $(fastp --version 2>&1 | sed -e "s/fastp //g")
END_VERSIONS
