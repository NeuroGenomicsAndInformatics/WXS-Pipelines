FROM ubuntu:22.04

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
	bzip2 \
	zlib1g \
	zlib1g-dev \
	openjdk-11-jdk \
	tabix \
	bcftools \
	rsync && \
	apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

ADD scripts /scripts

CMD ["/usr/bin/env bash"]
