
library(ggplot2)
load('figure3.related.Rdata')

## define a picture function

color.fun <- function(fn, var1, var2, colors, limit, bounder, na1, na2, panel){
  
  constant1 = grep(paste(sub('_rg','',var1),'R[1,2]_rg',sep=''),colnames(fn),value = T)
  constant2 = grep(paste(sub('_rg','',var2),'R[1,2]_rg',sep=''),colnames(fn),value = T)
  
  
  draw.fn <- fn[fn$n_sub<=limit & (!fn$geno %in% c('WT','non-functional')),]
  draw.fn = draw.fn[,c('geno',..var1,..var2,..constant1,..constant2)]
  ## color different quadrant
  draw.fn<- draw.fn %>%
    mutate(color = case_when(!!as.name(var1) > 1 & !!as.name(var2) > 1 ~ colors[1],
                             !!as.name(var1) < 1 & !!as.name(var2) > 1 ~ colors[2],
                             !!as.name(var1) < 1 & !!as.name(var2) < 1 ~ colors[3],
                             !!as.name(var1) > 1 & !!as.name(var2) < 1 ~ colors[4]))%>%
    ## variants that were not constantly in the same quadrant were colored as grey
    mutate(color = case_when(!!as.name(constant1[1]) > 1 & !!as.name(constant1[2]) < 1 ~ colors[5],
                             !!as.name(constant1[1]) < 1 & !!as.name(constant1[2]) > 1 ~ colors[5],
                             !!as.name(constant2[1]) > 1 & !!as.name(constant2[2]) < 1 ~ colors[5],
                             !!as.name(constant2[1]) < 1 & !!as.name(constant2[2]) > 1 ~ colors[5],
                             T ~ color))
  
  
  
  ## cal proportion of each color denoting four quadrant and those non-constant
  
  p.f = as.data.frame(t(table(draw.fn$color)/nrow(draw.fn)))
  
  p.f$Freq = round(p.f$Freq,3)
  
  p.f$Var2 = factor(p.f$Var2,levels = colors)
  
  p.f = p.f[order(p.f$Var2),]
  
  ## picturing
  
  p <- draw.fn%>%
    ggplot(aes(x=.data[[var1]],y=.data[[var2]]))+
    geom_point(aes(color= color),size=.6)+
    #geom_density2d(alpha=0.5, contour_var = "ndensity",color='black',show.legend = T)+
    scale_color_identity()+
    geom_vline(xintercept = 1,color='black',linetype='dashed')+
    geom_hline(yintercept = 1,color='black',linetype='dashed')+
    geom_text(x=1.4,y=1.5,label=paste(p.f$Freq[1]*100,'%',sep=''),color=p.f$Var2[1])+
    geom_text(x=-1.2,y=1.5,label=paste(p.f$Freq[2]*100,'%',sep=''),color=p.f$Var2[2])+
    geom_text(x=-1.3,y=-1.7,label=paste(p.f$Freq[3]*100,'%',sep=''),color=p.f$Var2[3])+
    geom_text(x=1.4,y=-1.7,label=paste(p.f$Freq[4]*100,'%',sep=''),color=p.f$Var2[4])+
    scale_x_log10(limits = c(-bounder,bounder),breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100),
                  oob = scales::squish)+
    scale_y_log10(limits = c(-bounder,bounder),breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100),
                  oob = scales::squish, position = 'right')+
    coord_equal(ratio = 1,expand = F)+
    theme_bw()+
    labs(x=expr(paste('Relative growth at ',!!na1,mu,"g/mL",sep="")),y=expr(paste('Relative growth at ',!!na2,mu,"g/mL",sep="")), title = panel)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(size = 10),
          axis.title = element_text(size=12),axis.text = element_text(size=10),title = element_text(size = 14))
  
  return(p)
  
}

## draw four panel using the above function

colors=c('#ff5f2e','#62578B','#478F62','#F0C05D','#a3a1a1')

p.A <- color.fun(fn, var1 = 'C0T3_rg',var2 = 'C2T3_rg',colors, limit = 1, bounder = 60, na1 = '0', na2 = '2',panel = 'A')

p.B <- color.fun(fn, var1 = 'C0T3_rg',var2 = 'C4T3_rg',colors, limit = 1, bounder = 60, na1 = '0', na2 = '4',panel = 'B')

p.C <- color.fun(fn, var1 = 'C0T3_rg',var2 = 'C2T3_rg',colors, limit = 6, bounder = 150, na1 = '0', na2 = '2',panel = 'C')

p.D <- color.fun(fn, var1 = 'C0T3_rg',var2 = 'C4T3_rg',colors, limit = 6, bounder = 150, na1 = '0', na2 = '4',panel = 'D')

library(patchwork)

p <- (p.A + p.B) / (p.C + p.D)


save(list = c('p','p.A','p.B','p.C','p.D','fn'), file = 'figure3.Rdata')
