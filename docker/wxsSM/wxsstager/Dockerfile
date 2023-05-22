##########################################
FROM ubuntu:22.04

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  bc \
  libz-dev \
	libbz2-dev \
  rsync \
  samtools \
  zlib1g \
  zlib1g-dev \
  liblzma-dev && \
  apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]
