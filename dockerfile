FROM ubuntu:20.04

RUN apt-get update
RUN apt install -y git
RUN cd /opt/
# this is private repo, make it public or the other way to copy the script into the container
RUN git clone https://github.com/TsailabBioinformatics/AGEseq2.git 
RUN python3 -m pip install Biopython pandas
RUN ["chmod", "+x", "/opt/AGEseq2/blat_binaries/blat_linux"]
RUN mkdir data

CMD python3 /opt/AGEseq2/AgeseqMain.py -t /data/targets.txt -r /data/reads/
