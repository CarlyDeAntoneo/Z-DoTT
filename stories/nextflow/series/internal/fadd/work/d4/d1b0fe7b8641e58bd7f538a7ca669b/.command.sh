#!/bin/bash -euo pipefail
cat input1/3281626_JCC510_Feb2025_KO_Input_IFN_B_S287_L004_R1_001.fastq.gz input3/3281626_JCC510_Feb2025_KO_Input_IFN_B_S287_L005_R1_001.fastq.gz input5/3281626_JCC510_Feb2025_KO_Input_IFN_B_S287_L006_R1_001.fastq.gz > KO_Input_IFN_B_1.merged.fastq.gz
cat input2/3281626_JCC510_Feb2025_KO_Input_IFN_B_S287_L004_R2_001.fastq.gz input4/3281626_JCC510_Feb2025_KO_Input_IFN_B_S287_L005_R2_001.fastq.gz input6/3281626_JCC510_Feb2025_KO_Input_IFN_B_S287_L006_R2_001.fastq.gz > KO_Input_IFN_B_2.merged.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:CAT_FASTQ":
    cat: $(echo $(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*$//')
END_VERSIONS
