#!/bin/bash

# DATE:2023-10-07
# NOTES: This is the code for Pacbio analysis
# Project: MCR-1
# AUTHOR: Lavenir-gzy
# Related Paper: An inevitable fitness cost facilitates the resistance reversibility in the absence of antibiotic selection(titled by 2023-10-07)


## Versions

# 1. samtools
    #samtools 1.6
    #Using htslib 1.6
# 2. blasr
    #blasr   5.3.5-5.3.5
# 3. ccs
    #ccs 6.4.0 (commit v6.4.0)
    #Using:
    #  unanimity : 6.4.0 (commit v6.4.0)
    #  pbbam     : 2.1.0 (commit v2.0.0-26-g05a8314)
    #  pbcopper  : 2.0.0 (commit v2.0.0-52-ga0c9454)
    #  boost     : 1.76
    #  htslib    : 1.15
    #  zlib      : 1.2.11
# 4. needle
    #EMBOSS:6.6.0.0

# 0. rawdata : m64061_210204_092329.subreads.bam

# 1. transfer rawdata to .sam file

conda activate pacbio # my conda environment built for pacbio analysis

samtools view -h -@ 50 -o m64061_210204_092329.subreads.sam m64061_210204_092329.subreads.bam

# 2. mapping to reference by BLASR to separate positive and negative subreads

  ## 2.1 blasr running

nohup blasr --bam --unaligned m64061_210204_092329.blasr.N1.unaligned.bam --out m64061_210204_092329.blasr.N1.bam m64061_210204_092329.subreads.bam ref_mcr_1.fasta> N1.blasr.log 2>&1 &
 
  ## 2.2 transfer blasr ouput into .sam format

#samtools view -h -@ 50 -o m64061_210204_092329.blasr.N1.sam m64061_210204_092329.blasr.N1.bam 

  ## 2.3 separate blasr ouput into two files(one for positive subreads and the other for negative subreads)

#sed -n '/^m/p' m64061_210204_092329.blasr.N1.sam | awk '$2 == 16 {print $1> "m64061_210204_092329.blasr.16.txt"} $2 == 0 {print $1> "m64061_210204_092329.blasr.0.txt"}'

samtools view -@ 50 ~/fitness_landscape/202205_New/Pacbio_analysis/m64061_210204_092329.blasr.N1.bam|awk '$2 == 16 {print $1> "m64061_210204_092329.blasr.16.txt"} $2 == 0 {print $1> "m64061_210204_092329.blasr.0.txt"} $2 == 256 {print $1> "m64061_210204_092329.blasr.256.txt"}'

 # our analysis uses primary alignment (FLAG = 16/0) only. And secondary alignment was ignored.

  ## 2.4 select raw subreads from rawdata(.sam)
time awk 'NR==FNR{a[$1]; next} $1 in a' m64061_210204_092329.blasr.16.txt m64061_210204_092329.subreads.sam >> m64061_210204_092329.n.subreads.sam

time awk 'NR==FNR{a[$1]; next} $1 in a' m64061_210204_092329.blasr.0.txt m64061_210204_092329.subreads.sam >> m64061_210204_092329.p.subreads.sam

  ## 2.5 transfer separated sam files into .bam format

samtools view -bS m64061_210204_092329.n.subreads.sam -o m64061_210204_092329.n.subreads.bam

samtools view -bS m64061_210204_092329.p.subreads.sam -o m64061_210204_092329.p.subreads.bam

  ## 2.6 add header

samtools view -H m64061_210204_092329.subreads.bam > header.sam

samtools reheader header.sam m64061_210204_092329.p.subreads.bam > m64061_210204_092329.p.reheader.subreads.bam

samtools reheader header.sam m64061_210204_092329.n.subreads.bam > m64061_210204_092329.n.reheader.subreads.bam

# 3. run ccs

  ## 3.1 ccs running

ccs --min-length 1400 --max-length 3200 -j 30 --min-passes 5 m64061_210204_092329.n.reheader.subreads.bam m64061_210204_092329.n.ccs.bam

ccs --min-length 1400 --max-length 3200 -j 30 --min-passes 5 m64061_210204_092329.p.reheader.subreads.bam m64061_210204_092329.p.ccs.bam

  ## 3.3 transfer ccs output into .sam format

samtools view -h -@ 50 -o m64061_210204_092329.p.ccs.sam m64061_210204_092329.p.ccs.bam

samtools view -h -@ 50 -o m64061_210204_092329.n.ccs.sam m64061_210204_092329.n.ccs.bam

# 4. call barcodes and genotype of ccs reads

python call.bar_geno.py

# 5. needle to get genotypes of ccs reads 

  ## 5.1 split call file into sub files for multi-processing

if [ ! -d ./split.call/ ]; then
echo "dir ./split.call/ doesn't exist"
mkdir -p ./split.call/
fi

awk 'BEGIN {n_seq=0;n_file=1} /^>/ {if(n_seq%80000==0){file=sprintf("./split.call/myseq%d.txt",n_file);n_file++} print >> file; n_seq++; next;} { print >> file; }' < m64061_210204_092329.call.bar_geno.txt

  ## 5.2 transfer .txt to .fa

if [ ! -d ./split.fa/ ]; then
echo "dir ./split.fa/ doesn't exist"
mkdir -p ./split.fa/
fi

cd ./split.call/
ls |while read i;do awk '/^>/ {printf("%s",">m64061_210204_092329/ccs/");sub(/[>]/,"",$0);sub(/[p]/,"/p",$0);sub(/[n]/,"/n",$0);printf("%s\n",$0);N++;next;} {sub(/[pn]_/,"",$1);printf("%s",$1);printf("\n")} END {printf("\n");}' $i > ./split.fa/${i::-4}.fa; done
cd ../

  ## 5.3 run needle
Rscript neele.r

  ## 5.4 genotypes with substitutions only were picked out
cd ./needle/
ls *.sam|while read i;do awk '/1626M/' $i | awk '{print $1,$4,$6,$10}' > ${i::-4}.call.blast.txt;done

  ## 5.5 get id-barcode relations from call files
cd ./split.call/
ls|while read i;do awk '/^>/ {printf("%s","m64061_210204_092329_ccs_");sub(/[>]/,"",$0);sub(/[p]/,"_p",$0);sub(/[n]/,"_n",$0);printf("%s ",$0);N++;next;} {sub(/[pn]/,"",$1);printf("%s",$2);printf("\n")} END {printf("\n");}' $i > ${i::-4}.call.barcode.txt;done

  ## 5.6 get mutations of reads by R
Rscript call.mut.r

cd ./mut/
ls | while read i;do cat $i >> ../merged.mut.txt;done

  ## file merged.mut.txt  includes cols: 1) reads ID; 2) results of needle(col 2;3); 3) sequence(col 4); 4) barcode(col 5) ; 5) mutations(col 6)

# 6. filter barcodes
  
  ## 6.1 only barcodes with right length were included: 15bp+15bp

awk '$5~/^[ATCG]{15}ATATGAGGCTTATCGTGAAG[ATCG]{15}$/' merged.mut.txt > merged.mut.15.txt

  ## 6.2 only barcode with no less than 2 reads were included

awk '{print $5}' merged.mut.15.txt |sort|uniq -c > uni_bar.txt

awk '$1 > 1' uni_bar.txt > filtered.uni_bar.txt

awk -F "\t" '{print $5”\t”$6}' merged.mut.15.txt |sort|uniq -c > bar_geno.txt

awk 'NR==FNR {a[$2]; next} $2 in a' filtered.uni_bar.txt bar_geno.txt > filtered.bar_geno.txt

  ## 6.3 barcodes with one substitution were regarded as the same 

    ## 6.3.1 split filtered.uni_bar.txt into sub-files for multi-processing 

awk 'BEGIN {n_seq=0;n_file=1} {if(n_seq%1900==0){file=sprintf("split.uni_bar%d",n_file);n_file++} print $2 >> file; n_seq++; next;} ' < ../filtered.uni_bar.txt

    ## 6.3.2 run Rscript to cluster similar barcodes(barcodes with no more than one substitutions)

Rscript simiar_barcode.py

    ## 6.3.3 combine outputs

find .| grep 'filtered'|while read i; do cat $i >> ../all.geno_bar.txt;done

awk -F “\t” '{print $1"\t"$2}' all.geno_bar.txt |sort|uniq > uniq.bar_geno.txt

awk '{gsub("ATATGAGGCTTATCGTGAAG","");print}' uniq.bar_geno.txt > uniq.bar_geno.30bp.txt

    ## 6.3.4 add barcodes of non-functional strains
awk '{print $1"\t"$2}' nofun.txt >> uniq.bar_geno.30bp.txt

    ## rename final output
mv uniq.bar_geno.30bp.txt final_bar_geno.txt

## END
## file final_bar_geno.txt is the final result of geno-barcode relations from pacbio sequencing
