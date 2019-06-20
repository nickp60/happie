FROM continuumio/miniconda3
RUN conda install -c bioconda cgview
RUN apt-get update
RUN apt-get install libgtk2.0-0:i386 libidn11:i386 libglu1-mesa:i386 -y
RUN which cgview
RUN which java
ENTRYPOINT [  "/opt/conda/bin/cgview" ]
