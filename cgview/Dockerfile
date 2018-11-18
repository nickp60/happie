FROM continuumio/miniconda3
RUN conda install -c bioconda cgview
RUN which cgview
RUN which java
ENTRYPOINT [  "/opt/conda/bin/cgview" ]
