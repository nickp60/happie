#! /bin/bash
set -e
source activate mobilephone


target=$1

function install_prophet {
    cd submodules/ProphET
    ./INSTALL.pl
    cd ../../
}

function update_prophet {
    echo "updating Prophet's prophage database"
    cd submodules/ProphET
    ./INSTALL.pl --update_db_only
    cd ../../
}

function test_prophet {
    # test prophet
    perl ./submodules/ProphET/ProphET_standalone.pl --fasta_in ./submodules/ProphET/test.fasta  --gff_in ./submodules/ProphET/test.gff --outdir t

}


function install_cafe {
    cd submodules/CAFE/
    make
    cd ../../
    ln -s submodules/CAFE/seg_clus seg_clus
}

function install_dimob {
	# # for dimob, we need a whole new env, and we need to fix some unescaped regexes with a horible sed call
	# conda create -n dimob perl-bioperl perl-log-log4perl perl-moose perl-config-simple hmmer
	# source activate dimob

	# #backup
	# cp ./submodules/islandpath/lib/Dimob/Config.pm ./submodules/islandpath/lib/Dimob/raw_Config.pm
    # cat ./submodules/islandpath/lib/Dimob/raw_Config.pm | sed s/'{{.+}}'/'\{\\\{.+\\\}\}'/g > ./submodules/islandpath/lib/Dimob/Config.pm
    echo "dimob no longer implemented"
    exit 1
}

function install_dependencies {
    # on osx legacy blast is not avaiable, you will need to get a compatible version from the biocore channel
    conda install -c biocore blast-legacy
    conda install perl-moose perl-test-requires perl-test-warnings perl-test-fatal perl-test-most perl-gd

    # deal with perl dependencies not available via conda

    cpanm install MooseX::Singleton --installdeps --force # The tests fail cause of a missing dependency.  There be dragons...
    cpanm install Bio::Coordinate
    cpanm install Bio::Graphics --force
    # cpanm install GD --force
    cpanm install GD::SVG --install-deps --force

}

function install_all {
    install_dependencies
    install_prophet
    install_cafe
    python setup.py develop

}

function main {
    case $1 in
	p)
	    install_prophet
	    ;;
	u)
	    update_prophet
	    ;;
	t)
	    test_prophet
	    ;;
	a)
	    install_all
	    ;;
	b)
	    install_dimob
            ;;
	c)
	    install_cafe
	    ;;
	d)
	    install_dependencies
	    ;;
	*)
	    echo "USAGE INSTALL.sh <a>"
	    ;;
    esac
}


main $target
