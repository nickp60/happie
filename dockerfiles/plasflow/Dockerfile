FROM python:3.5
MAINTAINER Nick Waters <nickp60@gmail.com>
RUN pip install https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.10.0rc0-cp35-cp35m-linux_x86_64.whl  numpy pandas scipy rpy2 scikit-learn biopython tzlocal
RUN git clone https://github.com/smaegol/PlasFlow
WORKDIR PlasFlow
RUN apt-get update && apt-get install r-base -y
RUN  R -e "source('http://bioconductor.org/biocLite.R'); biocLite('Biostrings')"
RUN pip uninstall  rpy2 -y
RUN pip install rpy2==2.8.3
# RUN apt-get install cpanminus -y
# RUN cpanm  Bio::Perl
#RUN cpanm GetOpt::longer
#ADD ecoli104.fna ./
#RUN ./PlasFlow.py  --input ecoli104.fna --output tsa
#RUN cat tsa
ENTRYPOINT [ "./PlasFlow.py" ]
