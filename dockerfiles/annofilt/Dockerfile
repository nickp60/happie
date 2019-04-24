FROM continuumio/miniconda3
MAINTAINER Nick Waters <nickp60@gmail.com>
RUN conda install -c bioconda blast
RUN git clone https://github.com/nickp60/annofilt.git ####
WORKDIR annofilt
RUN python setup.py install
RUN annofilt --version
ENTRYPOINT [ "annofilt" ]
