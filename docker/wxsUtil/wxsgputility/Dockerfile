##########################################
FROM nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  git \
  build-essential \
  zlib1g-dev \
  rsync \
  emacs \
  tabix \
  bedtools \
  python3 \
  openjdk-8-jre \
  && apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2.git \
  && cd bwa-mem2 \
  && make \
  && cp bwa-mem2* /usr/local/bin

RUN ln -s /usr/bin/python3 /usr/bin/python

CMD ["/bin/bash"]
