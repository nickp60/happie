FROM        ubuntu:14.04
MAINTAINER  Nick Waters nickp60@gmail.com
RUN apt-get update
RUN apt-get install \
        build-essential \
	ncbi-blast+ \
	default-jdk -y
RUN curl -o tmp.zip  -L https://downloads.sourceforge.net/project/brig/BRIG-0.95-dist.zip
RUN run unzip tmp.zip
