#!/bin/bash
set -euxo pipefail

FOLDER="$1"
PREMAP="$2"
BACKLINK="$(pwd)"

cd "$FOLDER"

# Required:
#  FOLDER:
#  - sequence.fa
#  - annotation.gff
#  PREMAP:
#  - pre-mapping.bed.gz
#  - pre-mapping.fa.gz

cp "${PREMAP}/pre-mapping.fa.gz" .
gunzip pre-mapping.fa.gz

cp "${PREMAP}/pre-mapping.bed.gz" .
gunzip pre-mapping.bed.gz

# Generate STAR index
STAR --runMode genomeGenerate --runThreadN "$(nproc)" \
  --genomeDir STAR-premap-index --genomeFastaFiles sequence.fa pre-mapping.fa \
  --sjdbGTFfile annotation.gff --sjdbOverhang 149 --limitGenomeGenerateRAM 128849018880 # 120GB RAM

cd "$BACKLINK"
