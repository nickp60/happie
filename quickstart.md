# Quickstart Guide
## Installation
To run happie, you will need to have Docker or Singularity installed, and python3.
 - [Install Docker](https://docs.docker.com/install/)
 - [Singularity](http://singularity.lbl.gov/)

Next, clone this repository, and install with python

```
git clone https://github.com/nickp60/happie
cd happie
python setup.py install
```

If you had any issues with the python install, try installing into a [bioconda](https://bioconda.github.io/) environment.

### Running happie the first time

Happie is now installed, but we need to get the docker/singularity images for all the programs it uses.  This is a bit of a pain, but having pre-built versions ensures that happie can be run replroducibly on any system, and can scale with distributed computing.

We have pre built all the images, and include a program, `happie_install_or_update` to do the work of downloading them from dockerhub.

Run that command:

```
happie_install_or_update
```

You should see the status as all the programs are installed.

## A simple analysis
### The data
Lets investigate the mobilomes of three /Staph/ strains:

- Staphylococcus aureus D30
- Staphylococcus aureus subsp. aureus MR1
- Staphylococcus aureus subsp. aureus str. CF-Marseille

First of all, lets make a directory, and download our genomes into it, and gunzip them.
```bash
mkdir staph
cd staph
curl -o D30.fasta.gz  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/171/455/GCF_000171455.1_ASM17145v1/GCF_000171455.1_ASM17145v1_genomic.fna.gz
curl -o MR1.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/176/195/GCF_000176195.1_ASM17619v1/GCF_000176195.1_ASM17619v1_genomic.fna.gz
curl -o Marseille.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/180/395/GCF_000180395.1_ASM18039v1/GCF_000180395.1_ASM18039v1_genomic.fna.gz
gunzip ./*
cd ../
```

You should have a directory with three genomes in it now.  You are ready to run happie.

### setting important paremeters
We know that staph genomes typically are between 2.4Mb and 3.4Mb, so we need to let happie know with `--QC_min_assembly` and `--QC_max_assembly`; it will throw an eerror for a genome outside the expected assembly size. For this analysis, lets skip running annofilt, the pipeline to remove bad annotations on the ends of contigs.  That requires a bit of work, so for now, we turn annofilt off wit `--skip_annofilt`.

We need to select one or more plasmid-detection tools.  Because mlplasmids doesnt have a built in database for Staph, lets use plasfow  by saying `--plasmid_tools plasflow`.

Lets say you're running on a machine with 4 cores; add `--cores 4` to the call

Finally, lets have it run abricate with vfdb to detect virulence factors, and the CARD database to detect AMR genes: `--analyses vfdb card`


We will run these commands in a bash loop, once for each genome:

```
for i in D30 MR1 Marseille;
do
   happie --contigs ./staph/${i}.fasta --QC_min_assembly 2400000 --QC_max_assembly 3400000 --skip_annofilt --plasmid_tools plasflow --analyses vfdb card --cores 4 --name ${i} --output happie_${i}
done
```

On my machine, this takes about 12 minutes per genome
