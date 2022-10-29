##########################################
FROM nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  git \
  build-essential \
  zlib1g-dev \
  rsync \
  python2 \
  python3 \
  && apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2.git \
  && cd bwa-mem2 \
  && make \
  && cp bwa-mem2* /usr/local/bin

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
