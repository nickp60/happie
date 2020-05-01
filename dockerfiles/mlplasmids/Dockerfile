FROM        r-base:latest
MAINTAINER  Nick Waters nickp60@gmail.com
RUN apt-get update
RUN apt-get install -y libgit2-dev build-essential libssl-dev libcurl4-openssl-dev libxml2-dev curl
RUN curl https://gitlab.com/sirarredondo/mlplasmids/raw/master/scripts/run_mlplasmids.R -o /bin/run_mlplasmids.R
RUN mkdir /data
RUN curl https://gitlab.com/sirarredondo/mlplasmids/raw/master/data/GCA_000250945.1_ASM25094v1_genomic.fna.gz -o /data/test.fasta.gz; gunzip /data/test.fasta.gz
RUN Rscript /bin/run_mlplasmids.R /data/test.fasta /data/test_results.tab .8
ENTRYPOINT ["/usr/bin/Rscript", "/bin/run_mlplasmids.R"]
