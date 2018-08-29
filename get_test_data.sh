#! /bin/bash

# to be run from the project root directory
mkdir test_data
cd test_data
for p in ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/299/455/GCA_000299455.1_ASM29945v1/GCA_000299455.1_ASM29945v1_genomic.fna.gz
do
    wget $p
done
gunzip ./*.gz
