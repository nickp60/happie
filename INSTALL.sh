#! /bin/bash

conda install -c biocore blast-legacy
cd submodules/ProphET
./INSTALL.pl
cd ../../
cpanm install MooseX::Singleton --force # The tests fail cause of a missing dependency.  There be dragons...
cpanm install Test::Most
cpanm install Bio::Coordinate
cpanm install Bio::Graphics
cpanm install GD
cpanm install GD::SVG

# for dimob, we need a whole new env, and we need to fix some unescaped regexes with a horible sed call
conda create -n dimob perl-bioperl perl-log-log4perl  perl-moose perl-config-simple hmmer
# backup
cp ./submodules/islandpath/lib/Dimob/Config.pm ./submodules/islandpath/lib/Dimob/raw_Config.pm
 cat ./submodules/islandpath/lib/Dimob/raw_Config.pm | sed s/'{{.+}}'/'\{\\\{.+\\\}\}'/g > ./submodules/islandpath/lib/Dimob/Config.pm
