#! /bin/bash
set -e
source activate mobilephone

case $1 in
    m)
	# on osx legacy blast is not avaiable, you will need to get a compatible version from the biocore channel
	conda install -c biocore blast-legacy
	conda install perl-moose perl-test-requires perl-test-warnings perl-test-fatal perl-test-most
	;;
    p)

	cpanm install MooseX::Singleton --installdeps --force # The tests fail cause of a missing dependency.  There be dragons...
	cpanm install Bio::Coordinate
	cpanm install Bio::Graphics
	cpanm install GD
	cpanm install GD::SVG

	cd submodules/ProphET
	./INSTALL.pl
	cd ../../
	;;
    u)
	cd submodules/ProphET
	./INSTALL.pl --updatedb
	cd ../../
	;;
    t)
	# test prophet
	perl ./submodules/ProphET/ProphET_standalone.pl --fasta_in ./submodules/ProphET/test.fasta  --gff_in ./submodules/ProphET/test.gff --outdir t
	;;
    h)

	python setup.py develop
	;;
    # d)

	# # for dimob, we need a whole new env, and we need to fix some unescaped regexes with a horible sed call
	# conda create -n dimob perl-bioperl perl-log-log4perl perl-moose perl-config-simple hmmer
	# source activate dimob

	# #backup
	# cp ./submodules/islandpath/lib/Dimob/Config.pm ./submodules/islandpath/lib/Dimob/raw_Config.pm
	# cat ./submodules/islandpath/lib/Dimob/raw_Config.pm | sed s/'{{.+}}'/'\{\\\{.+\\\}\}'/g > ./submodules/islandpath/lib/Dimob/Config.pm

    #       ;;

    c)
    # Install CAFE
	cd submodules/CAFE/
	make
	cd ../../
	;;
    *)
	echo "USAGE INSTALL.sh <a>"
	;;
esac
