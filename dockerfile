FROM ubuntu:20.04

RUN apt-get update
RUN apt install -y git
RUN apt-get -y install python3
RUN apt-get -y install pip
RUN python3 -m pip install Biopython 
RUN python3 -m pip install pandas
RUN apt-get -y install libcurl4
RUN git clone https://github.com/TsailabBioinformatics/AGEseq2.git
RUN cd AGEseq2; git checkout docker

CMD cd AGEseq2; python3 ./AgeseqMain.py -t /data/targets.txt -r /data/reads/
