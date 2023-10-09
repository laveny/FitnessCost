rm(list = ls())

setwd('/mnt/data3/disk/guoziyan/Pacbio/V2max/')

mcr1 = read.table('~/fitness_landscape/202205_New/Pacbio_analysis/ref_mcr_1.fasta')

mcr1 = as.character(mcr1[2,])

files <- list.files(path = './needle',pattern = '.call.blast.txt',full.names = T)

saving = '/mnt/data3/disk/guoziyan/Pacbio/V2max/mut/'

## 用于1对1比对的函数  这个函数估计有个bug 懒得改了,干脆不用了
#one_one <- function(seq,ref,start,end){
#  
#  ref = unlist(strsplit(ref,''))
#  
#  sequence = unlist(strsplit(sequence,''))[start:end]
#  
#  bwa = data.frame(ref,sequence)
#
#  bwa$position = c(1:1626)
#  
#  bwa = bwa[bwa$ref!= bwa$sequence,]
#  
#  mut = ''
#  
#  if (nrow(bwa) == 0) {
#    
#    mut = 'WT'
#    
#  } else {
#    
#    bwa$mut = paste(bwa$ref,bwa$position,bwa$sequence,sep = '')
#    
#    mut = paste(bwa$mut,collapse = ' ')
#  }
#  
#  return(mut)
#  
#  }
#
#


mclapply(files, function(x){
  
  fn <- data.table::fread(x, header = F)
  colnames(fn) <- c('id','start','map','seq')
  
  name = unlist(strsplit(x,'/'))
  name = name[length(name)]
  name = sub('.call.blast.txt','',name)
  
  ff <- data.table::fread(paste('/mnt/data3/disk/guoziyan/Pacbio/V2max/split.call/',name,'.call.barcode.txt',sep=''),header = F)
  colnames(ff) <- c('id','barcode')  
  fn <- merge(fn,ff,by = 'id')
  
  if ('TRUE' %in% is.na(fn)) {
    print(paste(name,' wrong 1'))
  }
  
  fn$mut =''
  fn = data.frame(fn)
  
  for (n in 1:nrow(fn)) {
    
    sequence = fn$seq[n]
    
    mapping = fn$map[n]
    
    id = fn$id[n]
    
    ## 因为事先挑出来了1626M的所有只有三种情况  2I1626M  1I1626M1I  1626M2I
    ## 另外写一个函数用于1对1比对
    ref = unlist(strsplit(mcr1,''))
    
    sequence = unlist(strsplit(sequence,''))
    
    if (mapping == '2I1626M') {
      
      sequence = sequence[3:1628]
      
    }else if (mapping == '1I1626M1I') {
      
      sequence = sequence[2:1627]
      
    }else if(mapping == '1626M2I'){
    
      sequence = sequence[1:1626]
      
    } else {print(paste(name,'wrong 2'))}
    
    
    
    bwa = data.frame(ref,sequence)
    
    bwa$position = c(1:1626)
    
    bwa = bwa[bwa$ref!= bwa$sequence,]
    
    mut = ''
    
    if (nrow(bwa) == 0) {
      
      mut = 'WT'
      
    } else {
      
      bwa$mut = paste(bwa$ref,bwa$position,bwa$sequence,sep = '')
      
      mut = paste(bwa$mut,collapse = ' ')
    }
    
    fn$mut[n] <- mut
    
  }
  
  write.table(fn,file = paste(saving,name,'.call.mut.txt',sep = ''),quote = F,sep = '\t',row.names = F,col.names = F)
  
},mc.cores = 31)
