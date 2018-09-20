#! /bin/bash
set -e

source activate mobilephone
# on osx legacy blast is not avaiable, you will need to get a compatible version from the biocore channel
conda install -c biocore blast-legacy
conda install perl-moose perl-test-requires perl-test-warnings perl-test-fatal perl-test-most

# cd submodules/ProphET
# ./INSTALL.pl
# cd ../../
cpanm install MooseX::Singleton --installdeps --force # The tests fail cause of a missing dependency.  There be dragons...
cpanm install Bio::Coordinate
cpanm install Bio::Graphics
cpanm install GD
cpanm install GD::SVG

# test prophet
perl ./submodules/ProphET/ProphET_standalone.pl --fasta_in ./submodules/ProphET/test.fasta  --gff_in ./submodules/ProphET/test.gff --outdir t


python setup.py develop
# for dimob, we need a whole new env, and we need to fix some unescaped regexes with a horible sed call
# conda create -n dimob perl-bioperl perl-log-log4perl perl-moose perl-config-simple hmmer

# source activate dimob

# backup
# cp ./submodules/islandpath/lib/Dimob/Config.pm ./submodules/islandpath/lib/Dimob/raw_Config.pm
#  cat ./submodules/islandpath/lib/Dimob/raw_Config.pm | sed s/'{{.+}}'/'\{\\\{.+\\\}\}'/g > ./submodules/islandpath/lib/Dimob/Config.pm
