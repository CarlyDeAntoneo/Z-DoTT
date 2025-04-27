#!/bin/bash -euo pipefail
printf "%s %s\n" KO_Z22_IFN_A_1.merged.fastq.gz KO_Z22_IFN_A_raw_1.gz KO_Z22_IFN_A_2.merged.fastq.gz KO_Z22_IFN_A_raw_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
    --quiet \
    --threads 2 \
    --memory 1024 \
    KO_Z22_IFN_A_raw_1.gz KO_Z22_IFN_A_raw_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_FASTP:FASTQC_RAW":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
