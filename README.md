[![DOI](https://zenodo.org/badge/146430799.svg)](https://zenodo.org/badge/latestdoi/146430799)
![icon](icon/logo.png)
# Happie
Happie is a wrapper for several independent programs used to identify different mobile elements. It is designed to thrive on HPC environments taking advantage of either Docker or Singularity to containerize the different programs.

Happie allows you to identify and extract all the "mobile" regions of a genome assembly -- those identified as being either prophage, plasmid, or genomic island sequence. This can be used to assess the charactaristics of the mobilome as opposed to the whole genome.

In the next few months, we will be uploaded a few case studies showing happie's utility.

(click here for a quickstart guide!)[quickstart.org]


## Install

Happie can manually be installed from github as follows

```
git clone https://github.com/nickp60/happie
cd happie
conda create -n happie  biopython pyyaml
source activate happie
python setup.py install
happie_install_or_update
```
That is going to take care of downloading Docker images for  [ProphET](https://github.com/jaumlrc/ProphET), [prokka](https://github.com/tseemann/prokka), [mlplasmids](https://gitlab.com/sirarredondo/mlplasmids), [dimob](https://www.brinkman.mbb.sfu.ca/~mlangill/islandpath_dimob/download.php), [mobsuite](https://github.com/phac-nml/mob-suite), [PlasFlow](https://github.com/smaegol/PlasFlow), and other goodies.  If you are running singularity, a folder called `.happie` will be created in your home directory.

## Running
(Optional) get some test data:

```
bash ./get_test_data.sh
```

And try running it!

```
happie --contigs ./test_data/ecoli104.fna  --output tmpresults
```


Or, if you are using singularity:

```
happie --contigs ./test_data/ecoli104.fna  --output tmpresults --virtualization singularity
```

You should get several output files:
```
tmp__small/
├── contig_names_key
├── happie_args.yaml
├── intermediate_files
│   ├── ProphET
│   ├── QC
│   ├── cgview.tab
│   ├── dimob
│   ├── mlplasmids
│   ├── mobile_abricate
│   ├── mobile_annofilt
│   ├── mobile_prokka
│   ├── wgs_abricate
│   ├── wgs_annofilt
│   └── wgs_prokka
├── results
│   ├── mobile_abricate.tab
│   ├── mobile_genome_coords
│   ├── mobile_test.fasta
│   ├── mobile_test.gbk
│   ├── mobile_test.gff
│   ├── wgs_abricate.tab
│   ├── wgs_test.gbk
│   └── wgs_test.gff
└── sublogs
    ├── QC.log
    ├── mobile_PropheET.log
    ├── mobile_abricate_resfinder_log.txt
    ├── mobile_abricate_vfdb_log.txt
    ├── mobile_annofilt.log
    ├── mobile_mlplasmids.log
    ├── mobile_prokka.log
    ├── wgs_abricate_resfinder_log.txt
    ├── wgs_abricate_vfdb_log.txt
    ├── wgs_annofilt.log
    └── wgs_prokka.log
```


##  Stages
Happie does the following things:
### 0) *QC*.
Happie removes short contigs and renames the contig names (to avoid issues with too-long names, non-standatd characters). The `names_key` file provides a link between the original names and the clean names.
### 1) Reannotates fasta with Prokka.
This ensures that the annotations are all the same format, which is useful for the pipelines that reference the annotations (like ProphET).

### 2) Identify Mobile Elements

#### Plasmids
If working with appropriate organisms, [[link][mlplasmids]] is your best. Otherwise, go with PlasFlow.  If you are feeling adventurous, try mob-suite, but use with caution.

#### Genomic Islands
Originally I was working on incorporating CAFE, but I later found out that it is not designed to handle certain organisms.  So, We went with IslandPath-DIMOB.

#### Prophages
See the short version of our head-to-head comparison here: [Testing 3 Prophage Finders](https://nickp60.github.io/weird_one_offs/testing_3_prophage_finders/).

###  3) Extract Mobile Elements

### 4) Assess features of Mobile Elements vs Chromosome
#### Abricate
#### VFDB
#### resfinder
#### AntiSmash

## FAQs
- Q: This seems like a lot of computational time could be saved by simply referencing the coordinates the mobile regions, rather than extracting and running analyses on the subset. A: you're right!
- Q: My harddrive is full!  why would you do this to me?  A: most of these tools rely on their own databases, which are all included in the docker images. Theres no way around it -- thats a loooooot of data.
- Q:

### Note
This module was renamed from "mobilephone", it was just too hard to google.


## Citing
After each run, happie creates a "citing.txt" file, with links to all the software happie uses.  Some of those use additional references, so make sure give credit to everyone involved!
