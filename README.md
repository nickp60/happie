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
We didn't install amos directly cause conda doesn't ahve it built for OSX :(.

```
conda create -n mobilephone emboss bedtools perl==5.22.0.1 perl-bioperl mob_suite perl-log-log4perl  perl-moose perl-config-simple hmmer
source activate mobilephone
# on osx legacy blast is not avaiable, you will need to get a compatible version from the biocore channel
conda install -c biocore bast-legacy
source deactivate; source activate mobilephone
./INSTALL.sh # this takes care of gettig the ProphET database,etc
```

cpanm install MooseX::Singleton
