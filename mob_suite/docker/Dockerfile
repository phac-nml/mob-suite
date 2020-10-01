FROM ubuntu:20.04
WORKDIR /
RUN apt update && apt -y upgrade && apt install -y python3-pip git curl libcurl4-openssl-dev libssl-dev
RUN git clone --single-branch -b dev https://github.com/phac-nml/mob-suite.git && cd mob-suite # && git checkout 3.0.0
RUN apt-get -y install ncbi-blast+ mash
RUN cd mob-suite && python3 setup.py install && rm -rf /mob-suite
RUN mob_init
CMD ["mob_typer", "-h"]
