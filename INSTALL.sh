#! /bin/bash

cd submodules/ProphET
./INSTALL.pl
cd ../../
cpanm install MooseX::Singleton --force # The tests fail cause of a missing dependency.  There be dragons...
