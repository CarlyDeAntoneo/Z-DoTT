#!/bin/bash -euo pipefail
cat input1/3281610_JCC510_Feb2025_WT_Input_Mock_A_S271_L004_R1_001.fastq.gz input3/3281610_JCC510_Feb2025_WT_Input_Mock_A_S271_L005_R1_001.fastq.gz input5/3281610_JCC510_Feb2025_WT_Input_Mock_A_S271_L006_R1_001.fastq.gz > WT_Input_Mock_A_1.merged.fastq.gz
cat input2/3281610_JCC510_Feb2025_WT_Input_Mock_A_S271_L004_R2_001.fastq.gz input4/3281610_JCC510_Feb2025_WT_Input_Mock_A_S271_L005_R2_001.fastq.gz input6/3281610_JCC510_Feb2025_WT_Input_Mock_A_S271_L006_R2_001.fastq.gz > WT_Input_Mock_A_2.merged.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:CAT_FASTQ":
    cat: $(echo $(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*$//')
END_VERSIONS
