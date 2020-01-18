FROM ubuntu:latest
MAINTAINER Nick Waters
# Version prokka:1.13b-2019
RUN apt-get update
RUN apt-get install curl locales build-essential libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl --yes
RUN export LC_ALL=en_US.UTF-8
RUN export LANG=en_US.UTF-8
RUN locale-gen en_US.UTF-8
RUN git clone https://github.com/tseemann/prokka.git
RUN curl  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz > ncbi.tar.gz && tar xzf  ncbi.tar.gz && mv ncbi-blast-2.10.0+/bin/* /usr/bin/
# get most recent :(
RUN curl ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz > tbl2asn.gz && gunzip tbl2asn.gz && chmod +x tbl2asn && mv tbl2asn /bin/ && rm prokka/binaries/linux/tbl2asn
RUN /usr/bin/perl /prokka/bin/prokka --setupdb
#RUN /usr/bin/perl /prokka/bin/prokka /prokka/test/plasmid.fna
ENTRYPOINT [ "/usr/bin/perl", "/prokka/bin/prokka" ]
