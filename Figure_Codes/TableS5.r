load('TableS5.related.Rdata')

var1 = 'C0T3_rg'
var2 = 'C4T3_rg'

constant1 = grep(paste(sub('_rg','',var1),'R[1,2]_rg',sep=''),colnames(fn),value = T)
constant2 = grep(paste(sub('_rg','',var2),'R[1,2]_rg',sep=''),colnames(fn),value = T)


fn <- fn[fn$n_sub<=6 & (!fn$geno %in% c('WT','non-functional')),]
fn = fn[,c('geno','AAchange','n_sub',..var1,..var2,..constant1,..constant2)]


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

rm(list = c('constant1','constant2','colors','confirm','label.fn','draw.fn'))

double.fn <- fn[fn$n_sub == 2 & fn$quadrant == 1  & fn$constant=='Y',]

double.fn = double.fn[,c('geno',..var1,..var2)]

double.fn <- double.fn%>%separate(geno,into = c('Gene1','Gene2'),remove = F)

double.fn$pos1 = as.integer(gsub('[ATCG]','',double.fn$Gene1))
double.fn$pos2 = as.integer(gsub('[ATCG]','',double.fn$Gene2))

single.fn <- fn[fn$n_sub == 1,]

## 1. from WT, in the situation with antibiotics, so conditions: 1 < WT < single mutant < double mutants  (C4T3_rg)

ref = as.numeric(wt.ref[,..var2][1])


result <- lapply(1:nrow(double.fn), function(x){
  
  
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
      
      if (index) {
        
        res = lapply(1:nrow(sing.wt), function(i){
          
          single = sing.wt$geno[i]
          sing.value = sing.wt[[var2]][i]
          
          df = data.frame(double_mutant = double, relative_growth_of_double = double.value, 
                          single_mutant = single, relative_growth_of_single = sing.value)
          return(df)
        })%>%rbind.fill()
        
        res$condition = 'WT;with antibiotics'
        return(res)
        
      }
      
    }  
    
  } 
  
  
}) %>% rbind.fill()



## 2）from single mutant, in the situation with antibiotics, so conditions: 1 < single mutant < double mutants (C4T3_rg)

result2.1 <- lapply(1:nrow(double.fn), function(x){
  
  
  double = double.fn$geno[x]
  double.value = double.fn[[var2]][x]
  
  sing1 = double.fn$Gene1[x]
  sing2 = double.fn$Gene2[x]
  
  sub = single.fn[single.fn$geno %in% c(sing1, sing2),]
  
  sing.wt = sub %>% filter(.data[[var2]] > (1 + 0.01))
  
  # whether rg of double mutant is > that of single mutant
  if (nrow(sing.wt) > 0) {
    
    index = double.value > (sing.wt[[var2]] + 0.01)
    
    if (T %in% index) {
      
      index = T
      
      if (index) {
        
        res = lapply(1:nrow(sing.wt), function(i){
          
          single = sing.wt$geno[i]
          sing.value = sing.wt[[var2]][i]
          
          df = data.frame(double_mutant = double, relative_growth_of_double = double.value, 
                          single_mutant = single, relative_growth_of_single = sing.value)
          return(df)
        })%>%rbind.fill()
        
        res$condition = 'single;with antibiotics'
        return(res)
        
      }
      
    }  
    
  } 
  
  
}) %>% rbind.fill()



## 3）from single mutant, in the situation without antibiotics, so conditions: 1 < single < double  (C0T3_rg)

result2.2 <- lapply(1:nrow(double.fn), function(x){
  
  
  double = double.fn$geno[x]
  double.value = double.fn[[var1]][x]
  
  sing1 = double.fn$Gene1[x]
  sing2 = double.fn$Gene2[x]
  
  sub = single.fn[single.fn$geno %in% c(sing1, sing2),]
  
  sing.wt = sub %>% filter(.data[[var1]] > (1 + 0.01))
  
  # whether rg of double mutant is > that of single mutant
  if (nrow(sing.wt) > 0) {
    
    index = double.value > sing.wt[[var1]] + 0.01
    
    if (T %in% index) {
      
      index = T
      
      if (index) {
        
        res = lapply(1:nrow(sing.wt), function(i){
          
          single = sing.wt$geno[i]
          sing.value = sing.wt[[var1]][i]
          
          df = data.frame(double_mutant = double, relative_growth_of_double = double.value, 
                          single_mutant = single, relative_growth_of_single = sing.value)
          return(df)
        })%>%rbind.fill()
        
        res$condition = 'single;without antibiotics'
        
        return(res)
        
      }
      
    }  
    
  } 
  
  
}) %>% rbind.fill()


result <- rbind(result, result2.1, result2.2)

all = ls()
all = all[!all %in% c('result')]
rm(list = c(all,'all'))

save.image('TableS5.Rdata')
