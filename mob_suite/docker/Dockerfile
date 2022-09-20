FROM ubuntu:20.04
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN apt update && apt install git python3-pip -y
RUN git clone https://github.com/phac-nml/mob-suite.git
RUN cd mob-suite #&& git checkout mob-3.0.3 && cd ..
RUN apt install libcurl4-openssl-dev libssl-dev -y
RUN pip3 install Cython numpy
RUN apt install mash ncbi-blast+ -y
RUN cd mob-suite && python3 setup.py install && cd .. && rm -rf mob-suite
RUN mob_init 
RUN apt clean
