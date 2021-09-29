FROM ubuntu:20.04

RUN apt-get update
RUN apt install -y git

RUN git clone https://github.com/TsailabBioinformatics/AGEseq2.git
RUN git checkout docker
RUN apt-get -y install python3
RUN apt-get -y install pip
RUN python3 -m pip install Biopython 
RUN python3 -m pip install pandas
RUN apt-get -y install libcurl4
CMD python3 ./AGEseq2/AgeseqMain.py -t /data/targets.txt -r /data/reads/
