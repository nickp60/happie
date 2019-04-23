FROM ubuntu:latest
MAINTAINER Nick Waters
RUN apt-get update
RUN apt-get install locales build-essential libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl --yes
RUN export LC_ALL=en_US.UTF-8
RUN export LANG=en_US.UTF-8
RUN locale-gen en_US.UTF-8
RUN git clone https://github.com/tseemann/prokka.git
RUN /usr/bin/perl /prokka/bin/prokka --setupdb
#RUN /usr/bin/perl /prokka/bin/prokka /prokka/test/plasmid.fna
ENTRYPOINT [ "/usr/bin/perl", "/prokka/bin/prokka" ]
