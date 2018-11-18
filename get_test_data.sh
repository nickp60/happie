#! /bin/bash

# to be run from the project root directory
mkdir test_data
cd test_data

curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/299/455/GCA_000299455.1_ASM29945v1/GCA_000299455.1_ASM29945v1_genomic.fna.gz -o ecoli104.fna.gz

curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/225/105/GCA_000225105.2_ASM22510v2/GCA_000225105.2_ASM22510v2_genomic.fna.gz -o ecoli_stec.fna.gz

gunzip ./*.gz
head -n 100 ecoli_stec.fna >  small.fasta
