#!/bin/bash

## DATE:2023-10-08
## This is the pipeline of ILLUMINA SEQUENCING data analysis
## AUTHOR: Lavenir-gzy

# 1. find barcodes from raw fastq reads

python barcode.calling.py

# 2. count the frequency of each barcode

python barcode.count.py

# 3. merge similar barcodes with one substitutiton (error)

  ## 3.1 use barcodes from pacbio analysis to create a pool of all possible similar barcodes
    
    ## 3.1.1 extract barcodes from pacbio analysis result
awk '{print $1}' final_bar_geno.txt > pacbio.all.bar.txt

    ## 3.1.2 split pacbio.all.bar.txt file into 100 sub-files for multi-processing
awk 'BEGIN {n_seq=0;n_file=1} {if(n_seq%1536==0){file=sprintf("split.uni_bar%d",n_file);n_file++} print $0 >> file; n_seq++; next;}' < pacbio.all.bar.txt

    ## 3.1.3 creating possible similar barcodes
Rscript create_all_sim_bar.run.r

    ## 3.1.4 combine outputs
find .|grep 'havesim'|while read i; do cat $i >> ../all_possible_sim_barcode.txt;done

    ## 3.1.5 get all barcodes from illumina sequencing and merge with that created from pacbio analysis result

python combine_all_illumina_barcode.py

awk 'NR==FNR {a[$0];next} $1 in a' "illumina.barcode.txt" "all_possible_sim_barcode.txt" > ../all_contained_possible.barcode.txt

# 4. calculate fitness by R

Rscript cal.fitness.R
