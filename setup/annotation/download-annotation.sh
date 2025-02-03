#!/bin/bash
set -euxo pipefail

ASSEMBLIES=$1

########################################################################################################################
# GRCm39
########################################################################################################################

# GENCODE
VERSION=M36

wget -O "${ASSEMBLIES}/GRCm39/GRCm39.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${VERSION}/GRCm39.primary_assembly.genome.fa.gz"

wget -O "${ASSEMBLIES}/GRCm39/gencode/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz"

wget -O "${ASSEMBLIES}/GRCm39/gencode/gencode.v${VERSION}.primary_assembly.annotation.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gtf.gz"

# RefSeq
wget -O "${ASSEMBLIES}/GRCm39/refseq/GRCm39.refseq.gff.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz"

########################################################################################################################
# CHM13 v2.0
########################################################################################################################

# Fasta
wget -O "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.gz" \
  "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"

# RefSeq
wget -O "${ASSEMBLIES}/CHM13v2/refseq/CHM13v2.refseq.gff.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"

# GRCh38 GENCODE annotation and sequence for lift-over
VERSION=47
wget -O "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/GRCh38.primary_assembly.genome.fa.gz"

wget -O "${ASSEMBLIES}/CHM13v2/gencode/liftoff/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz"

########################################################################################################################
# Postprocessing
########################################################################################################################

for file in ${ASSEMBLIES}/*/*.fa.gz
do
  gunzip "${file}"
  bgzip "${file/.fa.gz/.fa}"
  echo "Indexing" ${file}
  samtools faidx "${file}"
done
