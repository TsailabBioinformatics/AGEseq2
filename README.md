# Ageseq2
## Table of Contents

   * [Introduction](#introduction)
   * [Installation](#installation)
   * [Usage](#usage)
      * [Using the Command Line](#using-the-command-line)
      * [Using Docker](#using-docker)
   * [Information about Inputs](#information-about-inputs)
   * [Information about Outputs](#information-about-outputs)

### Introduction

Ageseq2-CLI is a python program for identification of target editing events from outputs of BLAT by aligning over thousands of reads to targets.

By default, the program will summarize the number of target editing events supported by the number of reads per target sequence in the provided file (see `Usage`).

Ageseq2 can be run on the command line interface or via Docker. This document will provide on instruction on how to use both. 

### Installation

#### Installation via Command Line
For advanced users, you need Python 3, Biopython, pandas, and BLAT binary to make Ageseq2 work.
You can find Python3 here: https://www.python.org/downloads/. Once python is installed in your computer, you can install both Biopython and pandas modules through pipy easily:

    python3 -m pip install Biopython
    python3 -m pip install pandas
 
In the package, BLAT binaries have been included in the blat_binaries folder, but cygwin1.dll might be required if you're running blat in a Windows system. You should be able to find it online with search engine. Once you found it, it has to stay with the other python scripts in the same place or the entry folder.
###### Possible issues with MAC
You might need to manually unlock the blat_macos in your system preference/setting.
###### Possible issues with LINUX
You might need to manually allow the blat_linux to be excuted with `chmod +x`.

Once the above is installed, clone this repository.
  
    git clone https://github.com/TsailabBioinformatics/AGEseq2/ 
    
- if cloning returned fatal: Authentication failed, then try again using a github personal access token. instructions to do so can be found here: https://docs.github.com/en/enterprise-server@3.1/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token
    
Once cloned, change directories into `AGEseq2`, and it should be ready to run. 

#### Installation via Docker
Install Docker <https://www.docker.com/products/docker-desktop>

Get the AGEseq image. Go to the docker image on a web browser:
<https://hub.docker.com/r/bendjamin101001/ageseq2>

Copy the 'Docker Pull Command', execute in the 'Terminal':
`docker pull bendjamin101001/ageseq2`

You should be able to see the 'bendjamin101001/ageseq' image in the docker dashboard, and it should be ready to run.

### Usage
The following documentation instructs how to run Ageseq2 on the command line or via Docker.

#### Using the Command Line  
1) Two inputs are required for Ageseq2:
    
    `python AgeseqMain.py -t [target_file] -r [reads_path] -sa [0|1]`

2) A target file and relative path to where reads files are stored are required for Ageseq. The default values are following:

    `target_file: targets.txt
    reads_path: ./reads`
    
3) Finally, you need to make sure the configuration file `AGEseq.conf` is in the same folder along with targets.txt. You should adjust those parameters to achieve a desired result.

4) If the target sequences are in `targets.txt`, and the reads files are in a folder named `reads`, and these two elements are in the same folder, then Ageseq can run with the simple command:

    `python AgeseqMain.py`
    
Additionally, `-sa` is set to 1 by default to not show alignments in the log file. If you're interested in looking each alignment of each read, you can change this to `-sa 0`.

###### Sapelo2 Specific
You will need to load pandas and Biopython with following commands before you can run Ageseq2. If Python3 is not automatically loaded, you can manually imported in a similar fashion:

    ml pandas/0.25.3-intel-2019b-Python-3.7.4
    ml Biopython/1.75-intel-2019b-Python-3.7.4
    
#### Using Docker
Open Docker and go to 'bendjamin101001/ageseq'

1) In the docker dashboard, click the blue ‘RUN▶️’ button.

2) Expand the optional setting.

3) Click the ‘…’ under Host Path, select the parent directory of ‘reads’ directory and ‘target.txt’

4) Enter ‘/data/’ as ‘Container Path’. You can enter a name for this run as ‘Container Name’.

5) Hit ‘Run’. Go to the ‘Containers / Apps’ tab on the left panel. Find the run you just launched. If successful, a CSV file containing a summary should be in the same directory as ‘reads’ and ‘target.txt’. If unsuccessful, check the ‘reads’ folder to make sure there are only .fastq files in there.

- If you did not specify a Container Name, it will be assigned as a random two words phrase. Don’t panic, just sort the containers by ‘Status’ or ‘Started Time’ to find it. The green logo means it is still running.

6) Click on the container to see the log. If you can see the progress running on the right side of the panel then everything is good. 

### Information about Inputs

##### Target format
A plain file with two columns, the first column is the name of target sequence, and the second column is the sequence. Now a fasta-style sequence file is also accepted as a valid target.
##### Reads format
Currently only fastq files are accepted. For example, merged amplicon sequences generated by PANDAseq in fastq can be used as reads files directly.
##### Configuration file
Paramters set by `AGEseq.conf`:

    remove_files      = 1   ;    # keep (0) or delete (1) intermediate files, default = 1
    READ_SNP_MINR_FREQ          = 0.05	;
    READ_INDEL_MINR_FREQ        = 0.001 ;
    READ_SNP_MINIMAL_SUPPORT    = 3 ;
    READ_INDEL_MINIMAL_SUPPORT  = 3 ;
    WOBBLE_BASE                 = True ; #Treat wobble base as one allele?
    WOBBLE_FREQ_LOW	            = 0.35 ; #Minimal frequency to call a wobble base
    WOBBLE_FREQ_HIGH            = 0.75 ; #Maximum frequency to call a wobble base

    #Below are paramters for BLAT

    blat_tileSize    = 7   ;    
    blat_oneOff      = 1   ;  
    blat_maxGap      = 20  ;    
    blat_minIdentity = 70  ;    
    blat_minScore    = 20  ; 


### Information about Outputs

In the summary table file, it includes following columns:

    input_file	targetID	aligned_target	aligned_consensus	sub_hits	editing_pattern_hits	editing_pattern

This summary file provides basic statistic information on targets and identified editing patterns from provide reads per file. Each identified editing pattern also has a consensus sequence built to show the identified location of gene edit.



