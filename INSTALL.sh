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
