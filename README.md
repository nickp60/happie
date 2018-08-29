# mobilephone
Mobilephone is a wrapper for several independent programs used to identify different mobile elements.

## Context independent Regions
These regions are interesting in and of themselves for what genes they carry, etc
### Phages
### Plasmids
### Genomic Islands
## Context-dependedt Regions
These regions are short, and are interestesitng solely because of their context
### Insertion Sequnes


# Install
For now, you have to make a conda env and then do a bit on manual stuff, cuase you have some databases to download.

```
git clone https://github.com/nickp60/mobilephone --recurse-submodules

for starters. the conda env:

```
git clone --recurse-submodules
conda create -n mobilephone
source activate mobilephone
conda install emboss bedtools
# on osx legacy blst is not avaiable, you will need to get a compatible version from the biocore channel
conda install -c biocore bast-legacy
source deactivate; source activate mobilephone
conda install perl-bioperl
git clone https://github.com/jaumlrc/ProphET; cd ProphET; ./INSTALL.pl

```
