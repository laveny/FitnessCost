

setwd('/mnt/data3/disk/guoziyan/MCR1/published/')

call_dir = './split.call/'

fa_dir = './split.fa/'

if (dir.exists(call_dir)) {
  
  print(paste("dir",call_dir,"has already exist",sep = ' '))
  
} else {
  
  if (dir.exists(fa_dir)) {
    
    print(paste("dir",fa_dir,"has already exist",sep = ' '))
    
  } else {
    
  dir.create(call_dir)
  dir.create(fa_dir)

  # 1. split call file into sub files for multi-processing
    
  system(command = "awk 'BEGIN {n_seq=0;n_file=1} /^>/ {if(n_seq%80000==0){file=sprintf(\"./split.call/myseq%d.txt\",n_file);n_file++} print >> file; n_seq++; next;} { print >> file; }' < m64061_210204_092329.max.call.txt")
  
  # 2. 
  
  system(command = "ls ./split.call/|while read i;do awk '/^>/ {printf(\"%s\",\">m64061_210204_092329/ccs/\");sub(/[>]/,\"\",$0);sub(/[p]/,\"/p\",$0);sub(/[n]/,\"/n\",$0);printf(\"%s\n\",$0);N++;next;} {sub(/[pn]_/,\"\",$1);printf(\"%s\",$1);printf(\"\n\")} END {printf(\"\n\");}' $i > ./split.fa/${i::-4}.fa; done")
  
  }
  
  
}

files <- list.files(pattern = 'fa$',path = fa_dir)

needle_dir = './needle'

reference = './ref_mcr_1.fasta'

mclapply(files, function(i){
  
  names = sub('.fa','',i)
  system(paste('cd ~/anaconda3/bin/;source ./activate;conda activate emboss;nohup time needle', reference, ' ',
               fa_dir,i,
               ' -gapopen 10 -gapextend 0.5 -aformat sam -outfile ',needle_dir,
               names,'.sam > ',needle_dir,names,'.log 2>&1 &',sep = ''))
  
},mc.cores = 31)

