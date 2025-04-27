#!/bin/bash -euo pipefail
cat input1/3281631_JCC510_Feb2025_KO_Z22_IFN_A_S292_L004_R1_001.fastq.gz input3/3281631_JCC510_Feb2025_KO_Z22_IFN_A_S292_L005_R1_001.fastq.gz input5/3281631_JCC510_Feb2025_KO_Z22_IFN_A_S292_L006_R1_001.fastq.gz > KO_Z22_IFN_A_1.merged.fastq.gz
cat input2/3281631_JCC510_Feb2025_KO_Z22_IFN_A_S292_L004_R2_001.fastq.gz input4/3281631_JCC510_Feb2025_KO_Z22_IFN_A_S292_L005_R2_001.fastq.gz input6/3281631_JCC510_Feb2025_KO_Z22_IFN_A_S292_L006_R2_001.fastq.gz > KO_Z22_IFN_A_2.merged.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:CAT_FASTQ":
    cat: $(echo $(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*$//')
END_VERSIONS
