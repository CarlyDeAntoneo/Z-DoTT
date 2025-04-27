#!/bin/bash -euo pipefail
printf "%s %s\n" WT_Input_Mock_C_1.merged.fastq.gz WT_Input_Mock_C_raw_1.gz WT_Input_Mock_C_2.merged.fastq.gz WT_Input_Mock_C_raw_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
    --quiet \
    --threads 2 \
    --memory 1024 \
    WT_Input_Mock_C_raw_1.gz WT_Input_Mock_C_raw_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTQC_RAW":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
