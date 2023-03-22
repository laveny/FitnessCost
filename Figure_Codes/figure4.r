
library(dplyr)


load('basefn.Rdata')

## to extract WT data

wt.fit = fn[fn$geno == 'WT',]


## only geno with no more than six substitutions were included

fn <- fn[fn$n_sub <= 6  & (! fn$geno %in% c('WT','non-functional')),]

## set vars

vars = paste('C',c('0','2','4'),'T3_fitness',sep = '')

## calculate fitness of each pre-mature codon

fn$Codon <- as.integer(fn$Codon)

pre.f<-fn %>% 
  group_by(Codon)%>%
  summarise(C0_mean = mean(!!as.name(vars[1])),C2_mean = mean(!!as.name(vars[2])),C4_mean = mean(!!as.name(vars[3])))


## set domains 

nums = c('D1','D2','D3') ##'transmembrane domain','linker','catalytic domain'
pre.f<-pre.f[pre.f$Codon<542 &  pre.f$Codon!=0,]
pre.f$domain = nums[2]
pre.f[pre.f$Codon<=175,'domain'] = nums[1]
pre.f[pre.f$Codon>=220,'domain'] = nums[3]



## plot figure

colors = c('#7E318E','#EC6200','#005E2A','#FFE200')

p1 <- pre.f %>%
  ggplot(aes(x=domain, y = C0_mean))+
  geom_jitter(alpha = 0.5,fill = 'grey',size = 1,shape = 16)+
  stat_boxplot(aes(color = as.character(domain)),geom ='errorbar', width = 0.4)+
  geom_boxplot(aes(group = domain, color = as.character(domain)),outlier.shape = NA,fill = NA)+
  geom_signif(test = 'wilcox.test',
              comparisons = list(c('D1','D2'),c('D2','D3'),
                                 c('D1','D3')),
              y_position = c(1.1,1.2,1.4),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0)+
  geom_point(aes(x='D4',y=as.numeric(wt.fit[,'C0T3_fitness'][1])),color='#D81159',shape=2,size=2)+
  scale_color_manual(values = c('D1'=colors[1],'D2'=colors[2],'D3'=colors[3]),
                     labels = c('transmembrane domain', 'linker', 'catalytic domain'))+
  scale_x_discrete(labels = c('transmembrane domain', 'linker', 'catalytic domain','wild-type'))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x='', y='Relative growth at \n0ug/mL colistin',title = 'A')+
  theme_classic()+
  coord_cartesian(ylim=c(0,1.5),clip = 'on')+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.title = element_text(size=12),axis.text = element_text(size=10),title = element_text(size = 14),
        legend.position = 'NA')

p2<-pre.f %>%
  ggplot(aes(x=domain, y = C2_mean))+
  geom_jitter(shape=16,alpha = 0.5,fill = 'grey',size = 1)+
  stat_boxplot(aes(color = as.character(domain)),geom ='errorbar', width = 0.4)+
  geom_boxplot(aes(group = domain, color = as.character(domain)),outlier.shape = NA,fill = NA,fill = NA)+
  geom_signif(test = 'wilcox.test',
              comparisons = list(c('D1','D2'),c('D2','D3'),
                                 c('D1','D3')),
              y_position = c(2.3,2.5,2.8),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0)+
  geom_point(aes(x='D4',y=as.numeric(wt.fit[,'C2T3_fitness'][1])),color='#D81159',shape=2,size=2)+
  scale_color_manual(values = c('D1'=colors[1],'D2'=colors[2],'D3'=colors[3]),
                     labels = c('transmembrane domain', 'linker', 'catalytic domain'))+
  scale_x_discrete(labels = c('transmembrane domain', 'linker', 'catalytic domain','wild-type'))+
  #scale_y_continuous(expand = c(0,0))+
  labs(x='', y='Relative growth at \n2ug/mL colistin',title = 'B')+
  theme_classic()+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.title = element_text(size=12),axis.text = element_text(size=10),title = element_text(size = 14),
        legend.position = 'NA')+
  coord_cartesian(ylim=c(0,4),clip = 'on')

p3<-pre.f %>%
  ggplot(aes(x=domain, y = C4_mean))+
  geom_jitter(shape=16,alpha = 0.5,fill = 'grey',size = 1)+
  stat_boxplot(aes(color = as.character(domain)),geom ='errorbar', width = 0.4)+
  geom_boxplot(aes(group = domain, color = as.character(domain)),outlier.shape = NA,fill = NA,fill = NA)+
  geom_signif(test = 'wilcox.test',
              comparisons = list(c('D1','D2'),c('D2','D3'),
                                 c('D1','D3')),
              y_position = c(2.3,2.5,2.8),map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),tip_length = 0
  )+
  geom_point(aes(x='D4',y=as.numeric(wt.fit[,'C4T3_fitness'][1])),color='#D81159',shape=2,size=2)+
  scale_color_manual(values = c('D1'=colors[1],'D2'=colors[2],'D3'=colors[3]),
                     labels = c('transmembrane domain', 'linker', 'catalytic domain'))+
  scale_x_discrete(labels = c('transmembrane domain', 'linker', 'catalytic domain','wild-type'))+
  labs(x='', y='Relative growth at \n4ug/mL colistin',title = 'C')+
  theme_classic()+
  theme(
    axis.title = element_text(size=12),axis.text = element_text(size=10),title = element_text(size = 14),
    legend.position = 'NA')+
  coord_cartesian(ylim=c(0,4),clip = 'on')


p <- p1 / p2 / p3

ggsave(filename = 'figure4.pdf',dpi=300,height = 16,width = 12,units = 'cm',plot = p)  


## compare with wild-type

wilcox.fn <- data.frame(Domain = rep(c('D1','D2','D3'),3),
                        WT = rep(c(wt.fit$C0T3_fitness,wt.fit$C2T3_fitness,wt.fit$C4T3_fitness),each = 3),
                        value = rep(c('C0_mean','C2_mean','C4_mean'),each = 3))

wilcox.fn$pvalue = 0
wilcox.fn$mean = 0
wilcox.fn$sig = ''

for (i in 1:9) {
  
  domain = wilcox.fn$Domain[i]
  
  var = wilcox.fn$value[i]
  
  mu = wilcox.fn$WT[i]
  
  X = unlist(pre.f[pre.f$domain == domain,var])
  
  Mean = mean(X)
  
  wilcox.fn$mean[i] = Mean
  
  a = wilcox.test(x = X,mu = mu)
  
  pvalue = a$p.value
  
  wilcox.fn$pvalue[i] <- pvalue
  
  if (pvalue < 0.001) {
    
    wilcox.fn$sig[i] = '***'
    
  } else if (pvalue < 0.01) {
    
    wilcox.fn$sig[i] = '**'
    
  } else if (pvalue < 0.05) {
    
    wilcox.fn$sig[i] = '*'
    
  } else {
    
    wilcox.fn$sig[i] = 'NS'
    
  }
  
}








save.image(file = 'figure4.Rdata') 



