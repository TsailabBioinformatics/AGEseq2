# Ageseq2
## Table of Contents

   * [Introduction](#introduction)
   * [Installation](#installation)
      * [Standard installation](#standard-installation)

   * [Usage](#usage)
      * [Simple Usage](#simple-usage)
      * [Inputs](#inputs)
      * [Outputs](#outputs)
   * [Citations](#citations)



### Introduction

Ageseq2-CLI is a command line interface (in Python) for identification of edited patterns from outputs of BLAT by aligning over thousands of reads to targets.

By default, the program will summarize the number of edited patterns supported by the number of reads per target sequence in the provided file (see `Usage`).


### Installation
#### Standard installation
For advanced users, you need Python 3.8+, Biopython, pandas, and BLAT binary to make Ageseq2 work.
You can find Python3 here: https://www.python.org/downloads/. Once python is installed in your computer, you can install both Biopython and pandas modules through pipy easily:

    python3 -m pip install Biopython
    python3 -m pip install pandas
 
In the package, BLAT binaries have been included in the blat_binaries folder, but cygwin1.dll will be required if you're running blat in a Windows system. You should be able to find it online with search engine. Once you found it, it has to stay with the other python scripts in the same place or the entry folder.
#### Possible issues with MAC
You might need to manually unlock the blat_macos in your system preference/setting.

### Usage

Two inputs are required for Ageseq2:
    
    python AgeseqMain.py -t [target_file] -r [reads_path] -sa [0|1]

A target file and relative path to where reads files are stored are required for Ageseq. The default values are following:

    target_file: targets.txt
    reads_path: ./reads
    
Finally, you need to make sure the configuration file `AGEseq.conf` is in the same folder along with targets.txt. You should adjust those parameters to achieve a desired result.

#### Simple Usage
This means that if you have put target sequences in the targets.txt, and reads files in a folder named "reads" in the same folder, you would be able to run Ageseq simply with  the following command:

    python AgeseqMain.py
    
Additionally, `-sa` is set to 1 by default to not show alignments in the log file. If you're interested in looking each alignment of each read, you can change this to `-sa 0`.

### Inputs
#### Target format
A plain file with two columns, the first column is the name of target sequence, and the second column is the sequence.
#### Configuration file
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


### Outputs
