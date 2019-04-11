FROM        ubuntu:14.04
MAINTAINER  Nick Waters nickp60@gmail.com
RUN apt-get update


# Installing prerequisites, Ubuntu 13.04:

### Use the Ubuntu Software Center or "sudo apt-get install" to install the following packages:
RUN apt-get install \
	build-essential \
	bioperl \
	git \
	libbatik-java \
	cpanminus \
	vcftools --yes
RUN git clone https://github.com/jonathancrabtree/Circleator
WORKDIR Circleator
RUN cpanm Log::Log4perl
RUN cpanm Text::CSV
RUN cpanm JSON
RUN cpanm SVG
RUN cpanm Test::Most
RUN cpanm Bio::FeatureIO::gff --force --verbose

RUN perl ./Build.PL
RUN ./Build
#RUN ./Build test
RUN ./Build install
RUN circleator --config=conf/genes-percentGC-GCskew-1.cfg --data=data/NC_011969.gbk >fig1.svg /install/path/bin/rasterize-svg fig1.svg png
ENTRYPOINT [ "/usr/local/bin/circleator" ]
