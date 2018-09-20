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

```
git clone https://github.com/nickp60/mobilephone --recurse-submodules
# or if already cloned
git submodule update --init
```

That is going to take care od downloading prophet as well.

For now, you have to make a conda env and then do a bit on manual stuff, cuase you have some databases to download.

After that, the conda env:

```
conda create -n mobilephone prokka  r-base r-r.utils r-devtools r-mlr emboss bedtools perl-extutils-pkgconfig  #perl==5.22.0.1 perl-bioperl mob_suite perl-log-log4perl  perl-moose perl-config-simple hmmer
./INSTALL.sh # this takes care of getting the ProphET database,etc
```