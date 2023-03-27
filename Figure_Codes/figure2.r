
load('figure2.related.Rdata')

WT.fit <- fn[fn$geno=='WT',]

fn = fn[fn$n_sub == 1,]

fn = fn[,c('geno','C0T3_rg','C2T3_rg','C4T3_rg')]

## 1. position of mutation

fn$pos = fn$geno

fn = separate(col = pos,into = c('ref','mut'),sep = '\\s*(\\d+)',convert = FALSE,data = fn)

fn$pos = fn$geno

fn$pos = sub('A','',fn$pos)
fn$pos = sub('T','',fn$pos)
fn$pos = sub('C','',fn$pos)
fn$pos = sub('G','',fn$pos)


## 2. domain

fn[fn$pos<=525,'domain'] = 'transmembrane domain'
fn[fn$pos>=660,'domain'] = 'catalytic domain'
fn[fn$pos>525 & fn$pos<660,'domain'] = 'linker'
unique(fn$domain)

### 3. function of significance 
fun_sig <-function(x){
  
  if (x<0.001) {
    
    y = '***'
    
  } else if (x<0.01) {
    
    y = '**'
    
  } else if (x<0.05) {
    
    y = '*'
    
  } else {
    y='ns'
  }
  
  
  return(y)
}


## 4.functions of picturing scatter(Panel A-C)

fun_fig.sca.sub <- function(data, pos, value,limit,color){
  
  ggplot(aes( x = .data[[pos]], y = .data[[value]]), data = data)+
    geom_point(aes(color = .data[[value]]), size = 0.2
    ) + 
    scale_color_gradientn(colours = color,limits = c(-limit,limit),oob=scales::squish,
                          guide = guide_colorbar(barwidth = 12,barheight = .6),
                          labels = as.character(10^c(-limit,0,limit)),breaks= c(-limit, 0 ,limit))+ 
    scale_x_continuous(breaks = seq(0,1626,300))+
    scale_y_continuous(
      breaks= c(-limit, 0 ,limit),limits = c(-limit,limit),labels = as.character(10^c(-limit,0,limit)),
      oob =scales::squish)+
    geom_hline(yintercept = 0,color='grey',linetype='dashed')+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
          axis.text = element_text(size=10),legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title.y = element_text(size=12,angle = 0, vjust = 0.5),legend.position = 'None',
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),
          title = element_text(size = 14),
          panel.background = element_blank(),plot.background = element_blank()
    )+
    labs(x='',y='')
  
}

fun_fig.sca <- function(fn,pos, value, color, limit){
  
  fn$pos = as.integer(fn$pos)
  
  fn[[value]] = log10(fn[[value]])
  
  pA <- fn %>% filter(mut == 'A')%>%
    fun_fig.sca.sub(pos,value,limit,color)+
    labs(x = 'Position')+
    theme(legend.position = 'top',axis.ticks.x = element_line(size = 1),axis.text.x = element_text(size = 10))+
    scale_x_continuous(breaks = seq(0,1626,300),position = 'top')
  
  
  pT <- fn %>% filter(mut == 'T')%>%
    fun_fig.sca.sub(pos,value,limit,color) 
  
  pC <- fn %>% filter(mut == 'C')%>%
    fun_fig.sca.sub(pos,value,limit,color) 
  
  pG <- fn %>% filter(mut == 'G')%>%
    fun_fig.sca.sub(pos,value,limit,color) 
  
  p <- pA / pT / pC / pG
  
  return(p)
  
}

# ----------------------------------------------------------------------------------------- #
#                                         Panel A-C                                         #
# ----------------------------------------------------------------------------------------- #

color_0 <- c('#053061','#3884bb','#a7d0e4','#f7f7f7','#f7b799','#ca4842','#67001f')

p.sca.0 <- fun_fig.sca(fn, pos = 'pos',value = 'C0T3_rg',color = color_0,limit = 1)

color_2 <- c('#276419','#6eae36','#c7e89f','#f7f7f7','#f5c4e1','#d6589e','#8e0152')

p.sca.2 <- fun_fig.sca(fn, pos = 'pos',value = 'C2T3_rg',color = color_2,limit = 2)

color_4 <- c('#2d004b','#715aa0','#bfbbda','#f7f7f7','#fdc57f','#d1740f','#7f3b08')

p.sca.4 <- fun_fig.sca(fn, pos = 'pos',value = 'C4T3_rg',color = color_4,limit = 2)


## 5.functions of calculating proportions(Panel D-F)


fun_cal.prop <- function(fn, vars , ref , cut.L, name.L){
  
  fn2 <- fn %>%
    dplyr::select(vars,domain,geno,pos) %>%
    ## grouping
    mutate( group = case_when(.[[vars]] < cut.L[1] ~ name.L[1],
                              .[[vars]] > cut.L[1] & .[[vars]] < cut.L[2] ~ name.L[2],
                              .[[vars]] > cut.L[1] ~ name.L[3]))
  
  ## group by domain and cal sum of mutations of each domain in each group
  ff <- fn2 %>% 
  group_by(domain, group) %>%
    summarise(count = n())%>%
    as.data.frame()%>%
    ## transfer sum into frequency
    group_by(domain)%>%
    mutate(prop := prop.table(count))
  
  
  
  return(list(fn2,ff))
  
}



## 6.functions of test(Panel D-F)

fun_bar.test = function(ff){
  
  ff <- data.frame(ff)
  
  a <- combn(unique(ff$domain),2)
  
  result <- lapply(1:3, function(x){
    
    d1 = as.character(a[1,x])
    
    d2 = as.character(a[2,x])
    
    f1 <- lapply(unique(ff$group),function(y){
      
      sub <- ff %>%
        filter(domain %in% c(d1,d2))%>%
        mutate(position = case_when( 
          group == y ~ 'isY',
          group != y ~ 'isN'))%>%
        as.data.frame()%>%
        group_by(domain, position)%>%
        summarise( num = sum(count))%>%
        as.data.frame()
      
      d1.isY =sub[sub$domain == d1 & sub$position=='isY','num']
      d1.isN =sub[sub$domain == d1 & sub$position=='isN','num']
      
      d2.isY =sub[sub$domain == d2 & sub$position=='isY','num']
      d2.isN =sub[sub$domain == d2 & sub$position=='isN','num']

      
      
      mytest2 <- min(phyper(d1.isY, d1.isY+d1.isN, d2.isY+d2.isN, d1.isY+d2.isY),
                     phyper(d2.isY, d2.isY+d2.isN, d1.isY+d1.isN, d1.isY+d2.isY))
      
      return(data.frame(group = y, domain1 = d1, domain2 = d2, d1.isY = d1.isY,
                        d1.isN = d1.isN, d2.isY = d2.isY, d2.isN = d2.isN, 
                        phyper.p = mytest2, sig.phyper = fun_sig(mytest2)))
      
    })%>%rbind.fill()
    
    return(f1)
    
  })%>%rbind.fill()
  
  return(result)
}

## 7.functions of picturing t(Panel D-F)
fun_bar.draw <- function(ff,colors){
  
  p <- ff%>%
    ggplot(aes(x = domain, y = prop))+
    geom_bar(aes(group = group, color = group),
             stat = 'identity',position = position_dodge(width = .6),
             width = 0.4,fill='white',size=.6)+
    scale_color_manual(values = colors)+
    theme_classic()+
    labs(x='',y='Proportion')+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.text = element_text(size=10),legend.title = element_blank(),
          axis.title = element_text(size=12),title = element_text(size = 14),
          legend.position = c(.92,.75),legend.direction = 'vertical',legend.key.size = unit(.2,'cm'),
          axis.title.x = element_blank(),legend.background = element_blank(),
          legend.text = element_text(size = 7))
  
  return(p)
}



## 8.functions of sliding windows & test(Panel G-I)


fun_Lchisq <- function(ss,windows){
  
  fn <- lapply(unique(ss$group), function(y){
    
    ff <- lapply(1:(1626-windows + 1), function(x){
      
      start = x
      end = start + windows -1
      
      sub <- ss %>%
        mutate(pos = as.integer(pos),
               if.pos = if_else(pos >= start & pos <= end , 'in.win','out.win'),
               if.group = if_else(group == y, 'in.group','out.group'))
      
      
      mydata <- table(sub$if.group,sub$if.pos)
      
      
      
      if ( ! 0 %in% mydata & nrow(mydata) == 2 & ncol(mydata) == 2) {
        
        test = chisq.test(mydata)
        
        pp = test$p.value
        
        FC <- (mydata['in.group','in.win']/sum(mydata[,'in.win'])) / (mydata['in.group','out.win']/sum(mydata[,'out.win']))
        
      }  else {
        
        pp = 1
        
        FC = 1
        
      }
      
      sig = fun_sig(pp)
      
      output = data.frame(group = y, position = mean(c(start, end)), FC = FC, chisq.p = pp,
                          sig = sig, rawdata = mydata)
      output <- unite(rawdata.Var1,rawdata.Var2,data = output,col = var)
      output <- dcast(output, group + position + FC + chisq.p + sig ~ var, value.var = 'rawdata.Freq')
      return(output)
      
    })%>%rbind.fill()
    
    ff$padjusted.bf <- p.adjust(ff$chisq.p,method = 'bonferroni')
    
    return(ff)
    
  })%>%rbind.fill()
  
  return(fn)
  
}


## 9.functions of picturing (Panel G-I)
fun_line.draw <- function(fn, colors){
  
  p.bf <- fn%>%
    mutate(pdraw = if_else(FC > 1, -log10(padjusted.bf),log10(padjusted.bf))) %>%
    ggplot(aes(x = position, y = pdraw, group = group))+
    geom_line(aes(color = group),size = .3)+
    scale_color_manual(values = c(colors), name = '')+
    theme_classic()+
    scale_x_continuous(breaks = seq(0,1626,300))+
    geom_hline(yintercept = c(log10(0.01),-log10(0.01)),linetype = 'dashed',
               color = 'grey')+
    scale_y_continuous(breaks = seq(-4,4,2),labels = c('-0.0001','-0.01','1','0.01','0.0001'))+
    theme(legend.direction = 'horizontal', legend.position = c(.4,.1),
          legend.box.background = element_blank(),
          axis.text = element_text(size = 10), axis.title = element_text(size =12),
          title = element_text(size = 14),legend.key.size = unit(.2,'cm'),
          legend.text = element_text(size = 7),legend.background = element_blank())+
    labs(y = 'Significance of \nenrichment', x = 'Position')
  
  
  p.line <- list(p.bf)
  names(p.line) <- c('p.bf')
  
  return(p.line)
  
} 



## 10.picturing (Panel G-I)
windows = 30
vars_0 = 'C0T3_rg'

ref_0 = WT.fit[[vars_0]]

cut.L_0 = c(ref_0,1)
name.L_0 = c('cost-increased', 'cost-reduced', 'costless')
prop_0 <- fun_cal.prop(fn, vars_0, ref_0, cut.L_0, name.L_0)
Lchisq.result_0 <- fun_Lchisq(prop_0[[1]] ,windows )
bar.test_0 <- fun_bar.test(prop_0[[2]])
prop_0[[2]]$domain <- factor(prop_0[[2]]$domain, levels = c('transmembrane domain','linker','catalytic domain'))
prop_0[[2]]$group <- factor(prop_0[[2]]$group, levels = c('cost-increased', 'cost-reduced', 'costless'))

p.bar_0 <- fun_bar.draw(prop_0[[2]],colors = c('#3884bb','#a7d0e4','#ca4842'))
p.bar_0 <- p.bar_0+
  geom_signif(y_position = c(0.55,0.6,0.65),xmin = c(0.8,1.8,0.8),xmax = c(1.8,2.8,2.8),
              annotation = c('***','***','*'),tip_length = 0,textsize = 3,size = .4,vjust = .5)+
  geom_signif(y_position = c(0.7,0.75,0.8),xmin = c(1,2,1),xmax = c(2,3,3),
              annotation = c('***','***','*'),tip_length = 0,textsize = 3,size = .4,vjust = .5)
p.line_0 <- fun_line.draw(Lchisq.result_0,colors = c('#3884bb','#a7d0e4','#ca4842'))  


p_0 <- p.sca.0/p.bar_0/p.line_0$p.bf+ plot_layout(heights = c(1,1,1,1,2,2))




vars_2 = 'C2T3_rg'

ref_2 = WT.fit[[vars_2]]
cut.L_2 = c(1,ref_2)
name.L_2 = c('non-resistant', 'resistance-decreased', 'resistance-enhanced')
prop_2 <- fun_cal.prop(fn, vars_2, ref_2, cut.L_2, name.L_2)
Lchisq.result_2 <- fun_Lchisq(prop_2[[1]] ,windows )
bar.test_2 <- fun_bar.test(prop_2[[2]])
bar.test_2[order(bar.test_2$group),]
prop_2[[2]]$domain <- factor(prop_2[[2]]$domain, levels = c('transmembrane domain','linker','catalytic domain'))
prop_2[[2]]$group <- factor(prop_2[[2]]$group, levels = c('non-resistant', 'resistance-decreased', 'resistance-enhanced'))
p.bar_2 <- fun_bar.draw(prop_2[[2]],colors = c('#6eae36','#f5c4e1','#d6589e'))
p.bar_2 <- p.bar_2+
  geom_signif(y_position = c(0.55),xmin = c(0.8),xmax = c(1.8),
              annotation = c('**'),tip_length = 0,textsize = 3,size = .4,vjust = .5)+
  geom_signif(y_position = c(0.65,0.7),xmin = c(1,2),xmax = c(2,3),
              annotation = c('***','***'),tip_length = 0,textsize = 3,size = .4,vjust = .5)+
  geom_signif(y_position = c(0.75,0.8),xmin = c(1.2,2.2),xmax = c(2.2,3.2),
              annotation = c('**','**'),tip_length = 0,textsize = 3,size = .4,vjust = .5)


Lchisq.result_2$group <- factor(Lchisq.result_2$group, levels = name.L_2)

p.line_2 <- fun_line.draw(Lchisq.result_2,colors = c('#6eae36','#f5c4e1','#d6589e'))  

p_2 <- p.sca.2/p.bar_2/p.line_2$p.bf+ plot_layout(heights = c(1,1,1,1,2,2))





vars_4 = 'C4T3_rg'

ref_4 = WT.fit[[vars_4]]
cut.L_4 = c(1,ref_4)
name.L_4 = c('non-resistant', 'resistance-decreased', 'resistance-enhanced')
prop_4 <- fun_cal.prop(fn, vars_4, ref_4, cut.L_4, name.L_4)
Lchisq.result_4 <- fun_Lchisq(prop_4[[1]] ,windows )
bar.test_4 <- fun_bar.test(prop_4[[2]])
bar.test_4[order(bar.test_4$group),]

prop_4[[2]]$domain <- factor(prop_4[[2]]$domain, levels = c('transmembrane domain','linker','catalytic domain'))
prop_4[[2]]$group <- factor(prop_4[[2]]$group, levels = c('non-resistant', 'resistance-decreased', 'resistance-enhanced'))
p.bar_4 <- fun_bar.draw(prop_4[[2]],colors = c('#715aa0','#fdc57f','#d1740f'))
p.bar_4 <- p.bar_4+
  geom_signif(y_position = c(0.7),xmin = c(2.2),xmax = c(3.2),
              annotation = c('*'),tip_length = 0,textsize = 3,size = .4,vjust = .5)

Lchisq.result_4$group <- factor(Lchisq.result_4$group, levels = name.L_4)

p.line_4 <- fun_line.draw(Lchisq.result_4,colors = c('#715aa0','#fdc57f','#d1740f'))  


p_4 <- p.sca.4/p.bar_4/p.line_4$p.bf+ plot_layout(heights = c(1,1,1,1,2,2))


p <- p_0 | p_2 | p_4


save.image('figure2.Rdata')
