FROM continuumio/miniconda3
MAINTAINER Nick Waters <nickp60@gmail.com>
RUN conda config --add channels defaults && \
	conda config --add channels conda-forge && \
	conda config --add channels bioconda
RUN conda install blast amos mash circlator mob_suite
# RUN git clone https://github.com/nickp60/mob0-suite ##
WORKDIR mob-suite
RUN ln -s conda-show-coords-path /usr/local/bin/show-coords
RUN mob_init
# ADD ecoli104.fna ./
# RUN mob_recon --infile ecoli104.fna --outdir my_out_dir
# RUN cat my_out_dir/contig_report.txt #
# RUN ls my_out_dir
# RUN head my_out_dir/chromosome.fasta
ENTRYPOINT [ "mob_recon" ]
