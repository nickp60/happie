# Standalone antiSMASH build
# VERSION 2018070601
# slightly modified be me to execute with external volumes
FROM antismash/standalone-lite:4.2.0
# the real maintainer
LABEL maintainer="Kai Blin <kblin@biosustain.dtu.dk>"
MAINTAINER Nick Waters <nickp60@gmail.com>
ENV ANTISMASH_VERSION="4.2.0"

# Grab the databases
WORKDIR /antismash-${ANTISMASH_VERSION}
RUN python download_databases.py

ADD instance.cfg antismash/config/instance.cfg

#VOLUME ["/input", "/output"]
#WORKDIR /output

ENTRYPOINT ["/usr/local/bin/run"]
