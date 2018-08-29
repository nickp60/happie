#! /bin/bash

cd submodules/ProphET
./INSTALL.pl
cd ../../

wget https://downloads.sourceforge.net/project/amos/amos-3.1.0-rc1.tar.gz
tar xzf amos*.tar.gz
cd amos-3.1.0-rc1
./configure
make
make install

pip install mob_suite -U
mob_init
