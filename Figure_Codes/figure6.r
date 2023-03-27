


## How the result was created

#F0 = seq(0,0.99,0.01)
#M = seq(0,0.99,0.01)
#
#cal1 = function(f0,m){
#  
#  p = 1 / 2
#  
#  ft = p * f0
#  
#  t = (log( 1 + (( 1 / p - 1) / (1 - f0 )))) / ( 1 - m)
#  
#  return(data.frame(M = m, F0 = f0, Ft = ft, T = t))
#  
#}
#
#
#fn = expand.grid(F0,M)
#
#result <- lapply(1:nrow(fn), function(i){
#  
#  return(cal1(fn$Var1[i],fn$Var2[i]))
#  
#})%>%rbind.fill(
  

## load result directly


load('figure6.related.Rdata')

# ----------------------------------------------------------------------------------------- #
#                            calculate real-world data                                      #
# ----------------------------------------------------------------------------------------- #

ft = c(36.25, 24.28, 1/44*100 ,3.7, 1/57*100, 1/20*100 , 14.15, 3.92, 3.66, 7.5, 5.03,
       2.56, 7.5, 8.33, 2.86, 1.11, 2.56, 6.67, 1/62 * 100, 0.57, 
       2.17, 13.43, 7.78, 8.17, 8.17, 8.57, 
       2.65, 2.74, 3.61, 8, 10.26, 6.1, 5.41)

ft <- as.numeric(ft)
Y = vector()
for (i in 1:length(ft)) {
  
  Y = append(Y, log(1/(ft[i] * 0.01) -1))
  
}


X = seq(from = 0, to = length(Y)-1) * 30 /45  ## the time in the real world is 45 times longer than that in the lab
result = lm(Y~X+0)
summary(result)

m.fit = (result$coefficients + log(1-ft[1] * 0.01)) / log(ft[1] * 0.01)

conf = confint(result)  

m.high = (conf[1] + log(1-ft[1] * 0.01)) / log(ft[1] * 0.01)
m.low = (conf[2] + log(1-ft[1] * 0.01)) / log(ft[1] * 0.01)

real.fit <- data.frame(mean = c(m.fit),
                       up = c(m.high),
                       low = c(m.low),
                       f0 = c(ft[1]/100))

rm(list = c('fit','i','m.fit','m.high','m.low','X','Y','conf','ft'))

# ----------------------------------------------------------------------------------------- #
#                            calculate reference lines                                      #
# ----------------------------------------------------------------------------------------- #

## for wt distribution


fn.barcode = fn.barcode[,c('barcode','geno','T0_readsSUM','C0T3_readsSUM')]

for (i in c('T0_readsSUM','C0T3_readsSUM')) {
  
  fn.barcode[,i] <- fn.barcode[,..i]/sum(fn.barcode[,..i]) 
  
}

WT.f <- fn.barcode[fn.barcode$geno=='WT',]

nofun.f <- fn.barcode[fn.barcode$geno=='non-functional',]%>%
  summarise_if(is.numeric,sum)

WT.f$C0T3_rg <- (WT.f$C0T3_readsSUM/WT.f$T0_readsSUM)/as.numeric(nofun.f$C0T3_readsSUM/nofun.f$T0_readsSUM)

rm(nofun.f)
rm(i)
rm(fn.barcode)

## for estimated relative growth of wt



wt.ref <- as.numeric(wt.ref$C0T3_rg)

## for variants


var.f = fn[fn$n_sub == 1,]

var.f <- var.f[order(var.f$C0T3_rg)]

var.0.05 <-quantile(var.f$C0T3_rg,0.05)

var.0.95 <-quantile(var.f$C0T3_rg,0.95)

var.mean = mean(var.f$C0T3_rg)

rm(fn)
# ----------------------------------------------------------------------------------------- #
#                                       picturing                                           #
# ----------------------------------------------------------------------------------------- #

colormap <- terrain.colors(10)

br = c(1,3,5,7,30,60,90,180,365,365*2,365*5)

br.real = br / 45  ## the time in the real world is 45 times longer than that in the lab


linetype.v <- data.frame(Y = c(var.0.05,var.0.95,var.mean,as.numeric(wt.ref)),Type = factor(c(2, 2, 1,4)),Color = c('#54534d','#54534d','#54534d','#EB472F'),stringsAsFactors = FALSE)


p3 <- WT.f%>%
  ggplot()+
  geom_histogram(aes(y = C0T3_rg), fill = '#FFC1BA')+
  theme_bw()+
  scale_x_continuous(expand = c(0,0), breaks = c(0,750,1500),limits = c(0,1500))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  labs(x = 'wild-type')+
  geom_hline(data = linetype.v %>% filter(Type == 4),aes(yintercept = Y, linetype = Type, colour = Color))+
  scale_linetype_manual(values = c(4),
                        labels = c('Estimated growth rate \nof wild-type'),
                        name='')+
  scale_color_identity()+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),panel.background = element_blank(),
        plot.background = element_blank(),legend.position ='NA',
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10))

p4 <- var.f%>%
  ggplot()+
  geom_histogram(aes(y = C0T3_rg),fill ='#BDBDBD')+
  theme_bw()+
  scale_x_continuous(expand = c(0,0),breaks =c(0,250,500),limits = c(0,500))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),oob=scales::squish)+
  labs(y = '', x = 'variants')+
  geom_hline(data = linetype.v %>% filter(Type != 4),aes(yintercept = Y, linetype = Type, colour = Color))+
  scale_linetype_manual(values = c(1,2),
                        labels = c('Average relative growth', '95 confidence interval'),
                        name='')+
  scale_color_identity()+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = 'NA',panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10)
  )





p5 <- ggplot(data = sim.f, aes(x = F0, y = M, fill = log10(T), z = log10(T)))+
  geom_tile(alpha = 0.8)+
  geom_contour(breaks = log10(br.real),color = 'white',linetype='solid')+
  scale_fill_gradientn(colours = colormap,breaks=log10(br.real),
                       labels = br, oob = scales::squish, name = 'Time(t)')+  #
  geom_point(data = real.fit, aes(x = f0 , y = mean))+
  geom_errorbar(data = real.fit, aes(x = f0, y=mean, ymin= low, ymax = up),width=.01)+
  geom_hline(data = linetype.v,aes(yintercept = Y, linetype = Type, colour = Color))+
  scale_linetype_manual(values = c(1,2,4),
                        labels = c('Average relative growth', '95 confidence interval','Estimated growth rate \nof wild-type'),
                        name='')+
  scale_color_identity()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x = 'Frequency of resistant strain (f0)',y='Relative growth of the resistant strain (m)', 
       title = 'B')+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 14),
        legend.position = 'NA',panel.background = element_blank(),
        plot.background = element_blank())+
  coord_equal()

pp <- p3 + p4 + plot_spacer() + plot_layout(widths = c(1,1,2)) 

pp <- p5 + pp

save.image('figure6.Rdata')


