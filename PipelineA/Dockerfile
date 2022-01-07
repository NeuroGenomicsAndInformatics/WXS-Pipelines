# This image is intended to be used as an end to end pipeline for a single sample from unmapped fastq to gvcf
# These values can be used to set version numbers later in the build, though will require resetting them

ARG ALPINE_VER=edge
ARG BWA_VER=0.7.17
ARG SAM_VER=1.9
ARG BED_VER=2.27.1
#ARG JAVA_VER=8-jre
##########################################
FROM alpine:${ALPINE_VER} AS builder

ARG BWA_VER
ARG SAM_VER
ARG BED_VER
ARG PICARD_VER

RUN apk add --no-cache \
	autoconf \
	automake \
	cmake \
	gcc \
	g++ \
	git \
	zlib-dev \
	bzip2-dev \
	libc-dev \
	libtool \
	make \
	musl \
	python2 \
	tini \
	wget \
	xz-dev \
	ncurses-dev

RUN git clone https://github.com/lh3/bwa.git \
	&& cd bwa \
	&& make

RUN wget -q -P /tmp https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 \
  && tar xjvf /tmp/samtools-${SAM_VER}.tar.bz2 -C /tmp \
  && cd /tmp/samtools-${SAM_VER}/ \
  && ./configure --without-curses \
  && make

RUN wget -q -P /tmp https://github.com/arq5x/bedtools2/releases/download/v${BED_VER}/bedtools-${BED_VER}.tar.gz \
  && tar xzvf /tmp/bedtools-${BED_VER}.tar.gz -C /tmp \
  && cd /tmp/bedtools${BED_VER%%.*}/ \
  && make

FROM alpine:edge AS java

RUN apk add --no-cache \
	zlib-dev \
	libstdc++ \
  openjdk8-jre \
	python2 \
  bash \
	&& mv /usr/lib/jvm/java-1.8-openjdk/jre/lib/amd64/jli/libjli.so /lib

ARG BWA_VER
ARG SAM_VER
ARG BED_VER

COPY --from=builder /bwa/bwa /usr/local/bin/
COPY --from=builder /tmp/samtools-${SAM_VER}/samtools /usr/local/bin/
COPY --from=builder /tmp/bedtools*/bin/bedtools /usr/local/bin/
COPY --from=builder /usr/lib/liblzma.so.5 /usr/lib/liblzma.so.5

COPY ./scripts/* /scripts/

#ENTRYPOINT ["/usr/bin/pipeA.bash"]
CMD ["/bin/bash"]
