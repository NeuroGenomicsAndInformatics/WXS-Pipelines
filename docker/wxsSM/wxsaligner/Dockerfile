##########################################
FROM nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1

RUN apt-get update && \
 	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
  rsync \
  python3 \
  && apt-get -y clean && \
	apt-get -y autoclean && \
	apt-get -y autoremove

COPY ./scripts/* /scripts/

CMD ["/bin/bash"]