FROM openjdk:11-jre

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
	bzip2 \
	zlib1g \
	zlib1g-dev \
	tabix \
	bcftools \
	rsync && \
	apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

ADD scripts /scripts

CMD ["/usr/bin/env bash"]
