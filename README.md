# Happie
Happie is a wrapper for several independent programs used to identify different mobile elements.

## Context independent Regions
These regions are interesting in and of themselves for what genes they carry, etc
### Phages
#### ProphET
### Plasmids
#### mlplasmids
### Genomic Islands
#### CAFE
## Context-dependedt Regions
These regions are short, and are interestesitng solely because of their context
### Insertion Sequnes
#### OASIS

# Install

```
git clone https://github.com/nickp60/happie --recurse-submodules
# or if already cloned
git submodule update --init
```

That is going to take care od downloading prophet as well.

For now, you have to make a conda env and then do a bit on manual stuff, cuase you have some databases to download.

After that, the conda env:

```
conda create -n happie prokka  r-base r-r.utils r-devtools r-mlr emboss bedtools perl-extutils-pkgconfig
./INSTALL.sh # this takes care of getting the ProphET database,etc
```

# Running



### Note
This module was renamed from "mobilephone", it was just too hard to google.
