FROM continuumio/miniconda3
RUN conda install -c bioconda cgview
RUN apt-get update -y
RUN apt-get install libgtk2.0.0 libidn11 libglu1-mesa -y
RUN which cgview
RUN which java
ENTRYPOINT [  "/opt/conda/bin/cgview" ]
