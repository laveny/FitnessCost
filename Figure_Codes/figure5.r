load('figure5.related.Rdata')

# ----------------------------------------------------------------------------------------- #
#                                              Panel A                                      #
# ----------------------------------------------------------------------------------------- #


var1 = 'C0T3_rg'
var2 = 'C4T3_rg'

constant1 = grep(paste(sub('_rg','',var1),'R[1,2]_rg',sep=''),colnames(fn),value = T)
constant2 = grep(paste(sub('_rg','',var2),'R[1,2]_rg',sep=''),colnames(fn),value = T)


fn <- fn[fn$n_sub<=6 & (!fn$geno %in% c('WT','non-functional')),]
fn = fn[,c('geno','AAchange','n_sub',..var1,..var2,..constant1,..constant2)]

## color different quadrant
fn<- fn %>%
  mutate(quadrant = case_when(!!as.name(var1) > 1 & !!as.name(var2) > 1 ~ 1,
                              !!as.name(var1) < 1 & !!as.name(var2) > 1 ~ 2,
                              !!as.name(var1) < 1 & !!as.name(var2) < 1 ~ 3,
                              !!as.name(var1) > 1 & !!as.name(var2) < 1 ~ 4))%>%
  ## variants that were not constantly in the same quadrant were colored as grey
  mutate(constant = case_when(!!as.name(constant1[1]) > 1 & !!as.name(constant1[2]) < 1 ~ 'N',
                              !!as.name(constant1[1]) < 1 & !!as.name(constant1[2]) > 1 ~ 'N',
                              !!as.name(constant2[1]) > 1 & !!as.name(constant2[2]) < 1 ~ 'N',
                              !!as.name(constant2[1]) < 1 & !!as.name(constant2[2]) > 1 ~ 'N',
                              T ~ 'Y'))
  
draw.fn <- fn %>%
  mutate(confirm = case_when(geno %in% confirm ~ 'Y',
                             T ~ 'N'))%>% 
  mutate(color = case_when(confirm == 'Y' ~ '#FF0018',
                           quadrant == 1 & confirm == 'N' ~ '#ff5f2e',
                           T ~ '#a3a1a1'))


  
  
## mark confirm genotype

label.fn = draw.fn[draw.fn$confirm == 'Y']
label.fn$AAchange = gsub(' ','+', label.fn$AAchange)

p.A  <- draw.fn%>%
  ggplot(aes(x=.data[[var1]],y=.data[[var2]]))+
  geom_point(aes(color= color,size = color,shape=confirm,alpha=color))+
  scale_color_identity()+
  scale_size_manual(values=c(.3,2,2))+
  scale_alpha_manual(values = c(.5,1,.5))+
  geom_vline(xintercept = 1,color='black',linetype='dashed')+
  geom_hline(yintercept = 1,color='black',linetype='dashed')+
  scale_y_log10(limits=c(0.3,130),
                breaks = c(0.3,0,1,10,100),labels = c(0.3,0,1,10,100),
                oob=scales::squish,position = "right")+
  scale_x_log10(limits=c(0.3,130),
                breaks = c(0.3,0,1,10,100),labels = c(0.3,0,1,10,100),
                oob=scales::squish)+
  ggrepel::geom_label_repel(aes(x = .data[[var1]], y = .data[[var2]], label = AAchange),data = label.fn,size = 2)+
  theme_bw(base_family = 'ArialMT')+
  labs(x=expression("Relative growth at 0"*mu*"g/mL"),
       y=expression("Relative growth at 4"*mu*"g/mL"),
       title = 'A')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = 'NA',
        axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 14))+
  coord_equal()


rm(list = c('constant1','constant2','colors','confirm','label.fn','draw.fn'))


# ----------------------------------------------------------------------------------------- #
#                                              Panel B                                      #
# ----------------------------------------------------------------------------------------- #

p.B <- expr.fn %>%
  ggplot(aes(x = factor(Mutant), y = MEAN, fill = time))+
  geom_errorbar( aes(x=factor(Mutant), ymin=MEAN-SD, ymax=MEAN+SD), position = position_dodge(width = 0.9), width=0.4, colour="#ADAEB3", alpha=0.6, size = .8)+
  geom_bar(position = position_dodge(width = 0.9), color = "black", stat = "summary", fun = "mean",width = 0.8)+
  geom_point(aes(y = value),position = position_dodge(width = 0.9),show.legend = F,size = 0.7)+
  scale_fill_manual(values = c('#00A499','#F29F7C'))+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  labs(x = '',y = 'Percentage of target \nstrains(%)',title = 'B')+
  theme( strip.text = element_blank(),
         strip.background = element_rect(color = 'NA',fill = 'NA',size = 10),strip.placement = 'outside',
         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
         axis.text = element_text(size = 10), axis.title = element_text(size = 12), 
         title = element_text(size = 14),
         legend.direction = 'horizontal',legend.position =c(.8,1.1),legend.key.size = unit(.1,'cm'),
         legend.title = element_blank(),legend.box.background = element_rect(colour = 'black'),
         plot.background = element_rect(fill = 'NA', colour = 'NA'))

# ----------------------------------------------------------------------------------------- #
#                                              Panel C                                      #
# ----------------------------------------------------------------------------------------- #
colnames(fn)



draw.fn <- fn %>%
  # only constant variants were included
  filter(constant == 'Y')%>%
  mutate(group = case_when(n_sub == 1 ~ 'single mutation',
                           n_sub >1 & n_sub <=6 ~'multiple mutations'))%>%
  group_by(group,quadrant)%>%
  summarise(count = n())

draw.fn <- draw.fn %>%
  group_by(group)%>% 
  mutate(SUM = sum(count))%>%
  filter(quadrant == 1)%>%
  data.frame()%>%
  dplyr::select(-'quadrant') %>%
  melt()

draw.fn$group <- factor(draw.fn$group, levels = c('single mutation','multiple mutations'))

p.C <- draw.fn%>%
  ggplot(aes( x = group, y = value))+
  geom_bar(aes(fill = variable),stat = 'identity',position = position_dodge(.7),width = .6)+
  scale_y_log10(expand = c(0,0),breaks = c(1,10,100,1000,10000))+
  theme_classic(base_family = 'ArialMT') +
  scale_fill_manual(values = c('SUM' = '#32373E','count' = '#F74831'), labels = c('all variants', 'costless resistant variants'),
                    name = '')+
  labs(y='Count', x = "Type of variants", title = 'C')+
  theme(legend.direction = 'horizontal',legend.position = c(.25,.95),
        legend.box.background = element_rect(color = 'black'),
        legend.text =  element_text(size = 10), axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        title = element_text(size = 14),legend.key.size = unit(.5, 'cm'))

# ----------------------------------------------------------------------------------------- #
#                                              Panel D.below                                      #
# ----------------------------------------------------------------------------------------- #

double.fn <- fn[fn$n_sub == 2 & fn$quadrant == 1  & fn$constant=='Y',]

double.fn = double.fn[,c('geno',..var1,..var2)]

double.fn <- double.fn%>%separate(geno,into = c('Gene1','Gene2'),remove = F)

double.fn$pos1 = as.integer(gsub('[ATCG]','',double.fn$Gene1))
double.fn$pos2 = as.integer(gsub('[ATCG]','',double.fn$Gene2))

p.D.below <- double.fn%>%
  ggplot(aes(x = 1:1626)) +
  scale_color_gradient2(midpoint=0,low = 'blue',high = 'red',mid = 'yellow',limits=c(-1,1),oob= scales::squish)+
  scale_fill_gradient2(midpoint=0,low = 'blue',high = 'red',mid = 'yellow',
                       oob= scales::squish,trans = 'log10', breaks = c(0,1,5,10,50,100),limits = c(1,150))+
  labs(fill = '', x = 'Position of mutation 1',
       y = 'Position of mutation 2')+
  geom_jitter(aes(pos1, pos2,fill = .data[[var2]]), shape=21,  size=2,alpha=.5,color='black') +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.text = element_text(size = 8),
        axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 0, "pt"))+
  xlim(1,1626)+
  ylim(1,1626)+
  coord_equal()

# ----------------------------------------------------------------------------------------- #
#                                              Panel D.above                                #
# ----------------------------------------------------------------------------------------- #

fun_sig <-function(x){
  
  if (x<0.001) {
    
    y = '***'
    
  } else if (x<0.01) {
    
    y = '**'
    
  } else if (x<0.05) {
    
    y = '*'
    
  } else {
    
    y='NS'
  }
  
  
  return(y)
}



fun_sliding <- function(win.fn, window, range ){
  
  pos.L = vector()
  fisher.L = vector()
  FC.L = vector()
  
  lapply(1:(range - window + 1), function(x){
    
    down_bound = x
    up_bound = x + window - 1
    
    win.fn$dis = 'out.win'
    win.fn[win.fn$pos1 >= down_bound & win.fn$pos1 <= up_bound, 'dis'] = 'in.win'
    win.fn[win.fn$pos2 >= down_bound & win.fn$pos2 <= up_bound, 'dis'] = 'in.win'
    
    mytable = table(win.fn$group, win.fn$dis)

    if ( ! (0 %in% mytable)){
      
      check = fisher.test(mytable)
      
      FC = (mytable['costless resistant','in.win']/sum(mytable['costless resistant',]))/(mytable['others','in.win']/sum(mytable['others',]))
      
      fisher.L <<- append(fisher.L, check$p.value)
      
      FC.L <<- append(FC.L, FC)
      
      pos.L <<- append(pos.L, mean(c(up_bound,down_bound)))
      
    } else {
      
      
      fisher.L <<- append(fisher.L, 1)
      
      FC.L <<- append(FC.L, 0)
      
      pos.L <<- append(pos.L, mean(c(up_bound,down_bound)))
      
    }
    
})
  
  result.fn <- data.frame( Position = pos.L,
                           fisher = fisher.L,
                           Fold_enrichment = FC.L)
  result.fn$fisher.sig = 0
  result.fn$chisq.sig = 0
  
  lapply(1:nrow(result.fn),function(x){
    
    result.fn$fisher.sig[x] <<- fun_sig(result.fn$fisher[x])
    result.fn$chisq.sig[x] <<- fun_sig(result.fn$chisq[x])
    
  })
  
return(result.fn)
  
}




sliding.f <- fn[fn$constant=='Y' & fn$n_sub == 2,] %>%
  mutate(group = case_when(quadrant == 1 ~ 'costless resistant',
                           T ~ 'others'))%>%
  separate(geno,into = c('Gene1','Gene2'),remove = F)

sliding.f$pos1 = as.integer(gsub('[ATCG]','',sliding.f$Gene1))

sliding.f$pos2 = as.integer(gsub('[ATCG]','',sliding.f$Gene2))

  
sliding.f <- fun_sliding(sliding.f, window = 50, range = 1626 )

p.D.above <- sliding.f%>%
  mutate(pdraw = if_else(Fold_enrichment > 1, -log10(fisher),log10(fisher))) %>%
  ggplot(aes(x = Position, y = pdraw))+
  geom_line(size = .3, aes(color = 'black'))+
  theme_classic()+
  scale_color_manual(name = "", values = c('fisher'='black'),labels = c("fisher's exact test"))+
  scale_x_continuous(breaks = seq(0,1626,300))+
  geom_hline(yintercept = c(log10(0.05),-log10(0.05)),linetype = 'dashed',
             color = 'grey')+
  scale_y_continuous(breaks = seq(-4,4,2),labels =math_format(expr = 10^-.x,format = abs),
                     limits = c(-4,4))+
  theme(legend.direction = 'horizontal', legend.position = c(.6,.8),
        legend.box.background = element_blank(),
        axis.text = element_text(size = 10), axis.title = element_text(size =12),
        title = element_text(size = 14),legend.key.size = unit(.2,'cm'),
        legend.text = element_text(size = 7),legend.background = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x = element_blank(),plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.title.x = element_blank()
  )+
  labs(y = 'Significance of ', x = '')

rm(list = c('expr.fn','draw.fn','sliding.f','fun_sig','fun_sliding'))

p.D <- plot_spacer() + p.D.above + plot_layout(heights = c(2,1)) 

p.D <- p.D / p.D.below


# ----------------------------------------------------------------------------------------- #
#                                              Panel E                                      #
# ----------------------------------------------------------------------------------------- #


single.fn <- fn[fn$n_sub == 1,]
 
## 1. from WT, in the situation with antibiotics, so conditions: 1 < WT < single mutant < double mutants  (C4T3_rg)

ref = as.numeric(wt.ref[,..var2][1])


result.fn <- lapply(1:nrow(double.fn), function(x){
  
  
  double = double.fn$geno[x]
  double.value = double.fn[[var2]][x]
  
  sing1 = double.fn$Gene1[x]
  sing2 = double.fn$Gene2[x]
  
  sub = single.fn[single.fn$geno %in% c(sing1, sing2),]
  
  # whether rg of sigle mutant is > wt
  sing.wt = sub %>% filter(.data[[var2]] > (ref + 0.01))
  
  # whether rg of double mutant is > that of single mutant
  if (nrow(sing.wt) > 0) {
    
    index = double.value > (sing.wt[[var2]] + 0.01)
    if (T %in% index) {
      
      index = T
      
    }  
    
  } else {
    
    index = F
    
  }
  
  ## index is the record of result
  if (index) {
    
    result = 'yes'
    
  } else {
    
    ## may be no or unknown
    
    if (nrow(sub) == 2) {
      
      result = 'no'
      
    } else {
      
      result = 'unknown'
    }
    
  }
  
  

  return(data.frame(gene = double, result = result))
  
}) %>% rbind.fill()

table(result.fn$result)

## result.fn
## no     yes unknown 
## 14       6       1 

## 2）from single mutant, in the situation with antibiotics, so conditions: 1 < single mutant < double mutants (C4T3_rg)

result2.1.fn <- lapply(1:nrow(double.fn), function(x){
  
  
  double = double.fn$geno[x]
  double.value = double.fn[[var2]][x]
  
  sing1 = double.fn$Gene1[x]
  sing2 = double.fn$Gene2[x]
  
  sub = single.fn[single.fn$geno %in% c(sing1, sing2),]
  
  sing.wt = sub %>% filter(.data[[var2]] > (1 + 0.01))
  
  if (nrow(sing.wt) > 0) {
    
    index = double.value > (sing.wt[[var2]] + 0.01)
    
    if (T %in% index) {
      
      index = T
      
    }  
    
  } else {
    
    index = F
    
  }
  
  if (index) {
    
    result = 'yes'
    
  } else {
    
    if (nrow(sub) == 2) {
      
      result = 'no'
      
    } else {
      
      result = 'unknown'
    }
    
  }
  
  
  
  return(data.frame(gene = double, result = result))
  
}) %>% rbind.fill()

table(result2.1.fn$result)

## result2.1.fn
## yes      no unknown 
##  11       9       1 

## 3）from single mutant, in the situation without antibiotics, so conditions: 1 < single < double  (C0T3_rg)


result2.2.fn <- lapply(1:nrow(double.fn), function(x){
  
  
  double = double.fn$geno[x]
  double.value = double.fn[[var1]][x]
  
  sing1 = double.fn$Gene1[x]
  sing2 = double.fn$Gene2[x]
  
  sub = single.fn[single.fn$geno %in% c(sing1, sing2),]
  
  sing.wt = sub %>% filter(.data[[var1]] > (1 + 0.01))
  
  if (nrow(sing.wt) > 0) {
    
    index = double.value > sing.wt[[var1]] + 0.01
    if (T %in% index) {
      
      index = T
      
    }  
    
  } else {
    
    index = F
    
  }
  
  ## index 是判断的结果
  if (index) {
    
    result = 'yes'
    
  } else {
    
   
    if (nrow(sub) == 2) {
      
      result = 'no'
      
    } else {
      
      result = 'unknown'
    }
    
  }
  
  

  return(data.frame(gene = double, result = result))
  
}) %>% rbind.fill()

table(result2.2.fn$result)


## result2.2.fn
## no unknown     yes
## 16       3       2

rm(list = c('sing1','sing2','ref','result','index','double','double.value','sub','single.wt','single.fn','wt.ref'))
rm(list = c('var1','var2','x'))

save.image('figure5.Rdata')
