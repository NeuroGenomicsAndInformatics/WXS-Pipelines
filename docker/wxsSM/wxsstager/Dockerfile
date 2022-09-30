FROM ubuntu:20.04

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
	unzip \
	bzip2 \
	zlib1g \
	zlib1g-dev \
	rsync && \
	apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

COPY ./scripts/* /scripts/

CMD ["/usr/bin/env bash"]
