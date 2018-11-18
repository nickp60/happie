# FROM continuumio/miniconda3
# RUN conda install -c bioconda abricate
# # test
# RUN abricate --check
# RUN abricate --list
FROM continuumio/miniconda3
MAINTAINER Nick Waters <nickp60@gmail.com>
RUN apt-get update
RUN apt-get install build-essential cpanminus emboss bioperl ncbi-blast+ gzip unzip \
	libjson-perl libtext-csv-perl libfile-slurp-perl liblwp-protocol-https-perl libwww-perl -y
RUN git clone https://github.com/tseemann/abricate.git
WORKDIR abricate
RUN cpanm install Path::Tiny
RUN ./bin/abricate --check
RUN ./bin/abricate --setupdb
RUN ./bin/abricate ./test/assembly.fa
