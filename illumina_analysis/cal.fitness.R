
## formal calculate fitness with two bio-repeats

## DATE : 2023.05.19


## 1. load in illumina sequencing data

rm(list = ls())

save_dir = './output/'

## located to where data saved

setwd('./Native/barcode.count/')

files <- list.files(pattern = '.txt')

## load in sample basis fn
basi <- data.table::fread("/mnt/data3/disk/guoziyan/MCR1/new20230516/sample_basis.csv")

## read barcode.count files

result <- lapply(files, function(x){
  
  id = sub('.count.txt','',x)
  
  name = basi[basi$sample_ID==id,]$official_ID
    
  fn <- read.table(x,header = F,sep='\t')
    
  colnames(fn) <- c('barcode',name)
    
  return(fn) 
  
}) 

## load in barcodes that should be remained
fn <- data.table::fread('./all_contained_possible.barcode.txt',header = F)

head(fn)

colnames(fn) <- c('barcode','main_barcode')

for (i in 1:length(result)) {
  
  if (result[i] != "NULL") {
    
    fn <- merge(fn, result[[i]], by = 'barcode',all.x = T)
    
  }
  
}

setwd(save_dir)

save(fn, file  = 'raw.reads.Rdata')  #2087611

## combine reads by main_barcodes

rm(list = ls())

load('raw.reads.Rdata')

fn[is.na(fn)]=0  ## 

fn <- fn %>% group_by(main_barcode) %>% select_if(is.numeric) %>%summarise_each(funs(sum)) 


## mapping to geno by pacbio analysis result
pacbio <- data.table::fread('final_bar_geno.txt',
                            header = F,stringsAsFactors = F,sep = "\t")

colnames(pacbio) = c('main_barcode','geno')

fn <- merge(fn, pacbio, by='main_barcode',all.x= T)


## output raw reads data 

setwd(save_dir)

save(fn, file  = 'raw.reads.geno.Rdata')


## sum up all bio-repeats

setDT(fn)

todo = unique(grep('T',gsub('R[123]','',colnames(fn)),value = T))

for (i in todo) {
  
  tosum = grep(paste(i,'R[123]',sep=''), colnames(fn), value = T)
  
  fn[,paste(i,'_readsSUM',sep='') := rowSums(.SD), .SDcols = tosum]
  
  
}


## 二、remove barcode has no reads at T0

fn <- fn[fn$`#T0_readsSUM` > 0,]



## 2. filtering WT barcodes

wt.f <- fn[fn$geno=='WT',]  



## condition 1: T0 sum > 100

cutoff = seq(5,200,5)

result <- lapply(cutoff, function(x){
  
  num <-wt.f%>%
    filter(`#T0_readsSUM`>=x)%>%
    nrow()
  num <- as.numeric(num)
  
  return(data.frame(cutoff = x, num = num))
  
})%>%rbind.fill()

p <-result %>% ggplot()+
  geom_point(aes(x = cutoff, y=num))+
  theme_bw()+
  labs(x = 'T0_readsSUM >= x', y = 'number of barcodes')

setwd(save_dir)
ggsave(filename = 'WT.barcode.minimum.reads.pdf',dpi=300,height = 12,width = 12,units = 'cm',plot = p,family='ArialMT') 

## set minimum reads cutoff = 100


wt.f <- wt.f[wt.f$`#T0_readsSUM` > 100,] 

## condition 2: MEAN +- SD

wt.f$C0T3.ratio = wt.f$C0T3_readsSUM/wt.f$`#T0_readsSUM`

MEAN = mean(wt.f$C0T3.ratio) 

SD = sd(wt.f$C0T3.ratio) 

sub0 <- wt.f %>% filter(C0T3.ratio >= MEAN-SD & C0T3.ratio <= MEAN + SD) 

wt.f$C2T3.ratio = wt.f$C2T3_readsSUM/wt.f$`#T0_readsSUM`

SD = sd(wt.f$C2T3.ratio) 

MEAN = mean(wt.f$C2T3.ratio)

sub2 <- wt.f %>% filter(C2T3.ratio >= MEAN-SD & C2T3.ratio <= MEAN + SD) 

wt.f$C4T3.ratio = wt.f$C4T3_readsSUM/wt.f$`#T0_readsSUM`

SD = sd(wt.f$C4T3.ratio) 

MEAN = mean(wt.f$C4T3.ratio)

sub4 <- wt.f %>% filter(C4T3.ratio >= MEAN-SD & C4T3.ratio <= MEAN + SD)

remain = intersect(intersect(sub0$main_barcode,sub2$main_barcode),sub4$main_barcode)

fn <- rbind(fn[fn$geno != 'WT',],fn[fn$geno == 'WT' & fn$main_barcode %in% remain,]) 

rm(sub0,sub2,sub4,wt.f,p,cutoff,i,MEAN,remain,SD,todo,tosum,result)


## filter variant
## for variants, sum of reads of all its barcode at T0 > 100
sub = fn[!fn$geno %in% c('WT','nofun'),] 

sub <- sub %>% dplyr::select(-main_barcode)%>%group_by(geno)%>%summarise_each(funs(sum))

cutoff = seq(5,200,5)

result <- lapply(cutoff, function(x){
  
  num <-sub%>%
    filter(`#T0_readsSUM`>=x)%>%
    nrow()
  num <- as.numeric(num)
  
  return(data.frame(cutoff = x, num = num))
  
})%>%rbind.fill()

p <- result %>% ggplot()+
  geom_point(aes(x = cutoff, y=num))+
  theme_bw()+
  labs(x = 'T0_readsSUM >= x', y = 'number of barcodes')

setwd(save_dir)
ggsave(filename = 'variant.minimum.reads.pdf',dpi=300,height = 12,width = 12,units = 'cm',plot = p,family='ArialMT') 

## set cutoff = 100

sub <-sub[sub$`#T0_readsSUM` > 100,]  

fn = rbind(fn[fn$geno %in% c('nofun','WT'),], fn[fn$geno %in% sub$geno,])  

rm(sub,p,result,cutoff)


## calculate fn_gene 和 fn_jitness

fn_gene <- fn %>% dplyr::select(-main_barcode)%>%group_by(geno)%>%summarise_each(funs(sum))

fn.fitness = fn_gene

fn.fitness <- as.data.frame(fn.fitness)

for (i in colnames(fn.fitness)) {
  
  if (!i %in% c('geno','#T0_readsSUM','#T0R1','#T0R2','#T0R3')) {
    
    fn.fitness <- fn.fitness %>% mutate(!!paste(sub('_readsSUM','',as.character(i)),'_rg',sep='') := fn.fitness[,i] /fn.fitness$`#T0_readsSUM`)
    
  }
  
}



for (i in grep('rg',colnames(fn.fitness),value = T)) {
  
  if (!i %in% c('geno','#T0_readsSUM')) {
    
    ref = as.numeric(fn.fitness[fn.fitness$geno=='nofun',i])
    
    fn.fitness[,i] = fn.fitness[,i]/ref
    
  }
  
}

setwd(save_dir)

save(list = c('fn_gene', 'fn.fitness'), file = 'basefn.Rdata')
