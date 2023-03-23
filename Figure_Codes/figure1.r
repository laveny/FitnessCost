library(ggplot2)
load('figure1.related.Rdata')

##  file introduction
## 1. fn was a file that each row is a genotype and columns containing 1) number of reads(_n_reads) of each group at different time point; 2) sum of reads of replicates(_readsSUM);
 3) estimated relative growth(_rg);4)n_subs means number of substitutions of that genotype; 5) AAchange denotes amino acid changes; 6) Codon is the position of the first stop codon
## 2. fn.barcode: each row denotes a barcode, whereas cols recording number of reads, similar to fn.
## 3. snr.fn: calculated results of snr. Col class denotes different concentration of Colistin. Col MEAN is the mean of SNR of each genotype within the group.
## 4. single mutant appeared in clinical isolates.

# ----------------------------------------------------------------------------------------- #
#                                             figure1. B                                    #
# ----------------------------------------------------------------------------------------- #

table(fn$n_sub)

textV = c('3091\n63.4%', 5914, 3158, 1208, 398, 109, 31)

color = c('#ff6433','#348c64','#a640a6','#fee342','#FD3469','#3E82C9','#757364')

count.f = data.frame(n_mutation = c(1:7),
                n_count = c(3091, 5914, 3158, 1208, 398, 109, 31),
                color = color)

p.B <- count.f%>% ggplot(aes(x=n_mutation,y=n_count,color = color))+
  geom_bar(stat = 'identity',width = 0.6,size=1,fill='white')+
  geom_text(aes(label=textV,color=color), vjust=-1,size=2.5)+
  scale_color_identity()+
  theme_classic()+
  labs(x='Number of substitutions',y='Number of variants',title = 'B')+
  scale_x_continuous(breaks = seq(1,7),labels = c('1','2','3','4','5','6','>6'))+
  scale_y_continuous(limits=c(0,8400),expand = c(0,0))+
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12),title = element_text(size = 14))

rm(list = c('textV','color'))

# ----------------------------------------------------------------------------------------- #
#                                             figure1. C                                    #
# ----------------------------------------------------------------------------------------- #


# transfer number of reads into frequency

a = grep('T',colnames(fn.barcode))

for (i in a) {
 
  fn.barcode[,i] = fn.barcode[,..i]/sum(fn.barcode[,..i]) 
   
}

a = grep('T3R[123]_n_reads',colnames(fn.barcode))

f.m <- fn.barcode[,..a]

f.m =cor(f.m,method = 'pearson')

p.C <- f.m %>% 
  as.data.frame() %>%
  rownames_to_column("f_id") %>%
  pivot_longer(-c(f_id), names_to = "samples", values_to = "counts") %>%
  ggplot(aes(x=samples, y=f_id, fill=counts)) + 
  geom_tile(color='white',size=.5) +
  scale_fill_gradientn(colours = brewer.pal(7,'YlOrRd'),limits=c(0.96,1),oob=scales::squish)+
  theme_get()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background =element_blank())+
  labs(x='Samples',y='Samples', title = 'C')+
  guides(fill= guide_colourbar(title = "R",barwidth = .7,barheight = 7))+
  #guides(fill=guide_legend(title='R'))+
  theme(legend.title = element_text(face='bold.italic'),legend.title.align=.8,
        axis.title = element_text(size=12),axis.text = element_text(size=10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),title = element_text(size = 14))+
  coord_equal()

for (i in c(2.5,4.5)) {
  p.C <- p.C+geom_segment(x=i,y=.5,xend=i,yend=6.5,color='black',size=.5)
}

for(i in c(2.5,4.5)) {
  p.C <- p.C+geom_segment(x=.5,y=i,xend=6.5,yend=i,color='black',size=.5)
}

rm(list = c('a','i'))

# ----------------------------------------------------------------------------------------- #
#                                             figure1. C.inset                              #
# ----------------------------------------------------------------------------------------- #

f.m['C0T3R1_n_reads','C0T3R2_n_reads']

p.C.inset <- fn.barcode %>%
  dplyr::select(C0T3R1_n_reads,C0T3R2_n_reads) %>%
  mutate_all(.funs = funs(.*1000000))%>%
   ggplot(aes(x=C0T3R1_n_reads,y = C0T3R2_n_reads))+
  geom_point(size=.005)+
  labs(x='replicate1',y='replicate2')+
  geom_text(aes(x = 1, y= 100000),label = "italic('R') ~ '>'~ '0.99'",parse = T)+
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x,n = 3),
    labels = trans_format("log10", math_format(10^.x)))+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x,n = 3),
    labels = trans_format("log10", math_format(10^.x)))+
  coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size=4),axis.text = element_text(size=3),title = element_text(size = 5))+
  geom_abline(slope = 1,color='grey',linetype='dashed',size=.4)
  

# ----------------------------------------------------------------------------------------- #
#                                             figure1. D                                    #
# ----------------------------------------------------------------------------------------- #


p.D <- snr.fn%>%ggplot(aes(x=group,y=MEAN,group=class,fill = class))+
  geom_errorbar(aes(ymin = MEAN-SE,ymax = MEAN+SE,group=class),position = position_dodge(.8),
                stat="identity",width=0.5,show.legend = F,color = 'black')+
  geom_bar(stat = 'identity',position = position_dodge(.8),width = 0.8,color = 'black')+
  scale_fill_manual(values = c('#DE3735','#FEA316','#58662B'),#c('#E79C4C','#214A5B','#7F3235'),,'#FDBE83','#C8A3B5','#2F4E58'
                    name=NULL,labels=c(expression("0"~""*mu*"g/mL"),
                                       expression("2"~""*mu*"g/mL"),expression("4"~""*mu*"g/mL")))+
  scale_y_log10(breaks = c(1,2,5,10,20,40), expand=c(0,0))+
  labs(x='',y='Signal-to-noise ratio', title = 'D')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size=12),axis.text = element_text(size=10),
        title = element_text(size = 14),
        legend.direction = 'horizontal',legend.position = c(.75,.9))
  
# ----------------------------------------------------------------------------------------- #
#                                             figure1. E                                    #
# ----------------------------------------------------------------------------------------- #


sig.fun <- function(p.value){
  
  if (p.value<0.001) {
    pv = '***'
  } else if (p.value<0.01) {
    pv = '**'
  } else if (p.value<0.05) {
    pv = '*'
  } else {
    pv = 'NS'
  }
  return(pv)
}

cal.chisq.fun <- function(fn, cutoff, value, fset){
  
  fn <- fn %>% 
    mutate(window = case_when(geno %in% fset ~ 'in',
                                   T ~ 'out'),
           value = case_when(!!as.name(value) > cutoff ~ 'True',
                             T ~ 'False'))
  
  mydata = table(fn$window,fn$value)
  
  test = chisq.test(mydata)
  
  sig = sig.fun(test$p.value)
  
  FC = (mydata['in','True']/sum(mydata['in',]))/(mydata['out','True']/sum(mydata['out',]))
  
  return(data.frame(cut = cutoff, p.value = test$p.value, sig = sig, Fold_enrichment = FC))
  
}

## only variants with no more than six substitutions were included

cli.f <- fn %>% filter(n_sub <=6 & (! geno %in% c('WT','non-functional')))

## grep variants that contain single variants in clinical-related isolates

cli.sin = unique(unlist(lapply(cli.sin,function(i) grep(i,fn$geno,value = T))))

cli.f <- lapply(seq(2,18,3), function(x){
  
  ff <- cal.chisq.fun(cli.f, x , 'C4T3_rg', cli.sin)
  
  return(ff)
  
}) %>% rbind.fill()

colors = c('#FFCD35','#FFCD35','#F7970E','#FC7B18','#F26600','#E95700')
cli.f$colors = colors[1:nrow(cli.f)]

p.E <- cli.f%>%
  ggplot(aes(x=cut,y=Fold_enrichment,fill = colors))+
  geom_bar(stat = 'identity')+
  scale_fill_identity()+
  geom_hline(yintercept = 1,linetype='dashed', color = '#898C96')+
  geom_text(aes(label= sig),size=4,vjust=-.1)+
  theme_bw()+
  labs(x= expression(paste('Threshold for relative growth(',phantom(0) >= "x)",sep='')),y='Fold enrichment for \nvariants in clinical isolates')+
  scale_x_continuous(breaks = c(seq(2,18,3)))+
  scale_y_continuous(limits = c(0,2.3))+
  labs(title = 'E')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_text(size=10),axis.title = element_text(size=12),title = element_text(size = 14))

rm(list = c('cli.sin','colors','cal.chisq.fun','sig.fun'))

save.image('figure1.Rdata')

