# Standalone antiSMASH build
# VERSION 20190618
# slightly modified be me to execute with external volumes
FROM antismash/base5:5.0.0rc1
# the real maintainer
LABEL maintainer="Kai Blin <kblin@biosustain.dtu.dk>"
# the maintainer of this  particular file
MAINTAINER Nick Waters <nickp60@gmail.com>
ENV ANTISMASH_VERSION="5.0.0"

ENV ANTISMASH_URL="https://dl.secondarymetabolites.org/releases/"
ENV ANTISMASH_VERSION="5.0.0"

ENV LANG C.UTF-8

# Grab antiSMASH
RUN curl -L ${ANTISMASH_URL}/${ANTISMASH_VERSION}/antismash-${ANTISMASH_VERSION}.tar.gz > /tmp/antismash-${ANTISMASH_VERSION}.tar.gz && \
    tar xf /tmp/antismash-${ANTISMASH_VERSION}.tar.gz && \
    rm /tmp/antismash-${ANTISMASH_VERSION}.tar.gz

ADD instance.cfg /antismash-${ANTISMASH_VERSION}/antismash/config/instance.cfg


RUN pip3 install /antismash-${ANTISMASH_VERSION} && python3 -c "import antismash; antismash.config.build_config(['--databases', 'mounted_at_runtime'], modules=antismash.get_all_modules()); antismash.main.prepare_module_data()"

WORKDIR /usr/local/bin
RUN ln -s /antismash-${ANTISMASH_VERSION}/docker/run

RUN mkdir /matplotlib && MPLCONFIGDIR=/matplotlib python3 -c "import matplotlib.pyplot as plt" && chmod -R a+rw /matplotlib

WORKDIR /

ADD instance.cfg /usr/local/lib/python3.5/dist-packages/antismash/config/instance.cfg
RUN mkdir -p /databases && download-antismash-databases --data /databases

# RUN find / -name  run_antismash.py
ADD small.fasta ./
RUN /antismash-5.0.0/run_antismash.py small.fasta --genefinding-tool prodigal -v --output-dir tmpoutput && rm -r tmpoutput
ENTRYPOINT [ "/antismash-5.0.0/run_antismash.py" ]
