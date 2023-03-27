
load('TableS2_4.related.Rdata')

fn = fn[fn$n_sub == 1 | fn$geno == 'WT',]

sub1 = fn[,c('geno','C0T3_rg')]

setwd('~/fitness_landscape/202205_New/Native_Promoter/github/upload/')

write.table(sub1, file = 'TableS2.txt', quote = F, row.names = F, col.names = F, sep = '\t')


sub2 = fn[,c('geno','C2T3_rg')]

setwd('~/fitness_landscape/202205_New/Native_Promoter/github/upload/')

write.table(sub2, file = 'TableS3.txt', quote = F, row.names = F, col.names = F, sep = '\t')


sub3 = fn[,c('geno','C4T3_rg')]

setwd('~/fitness_landscape/202205_New/Native_Promoter/github/upload/')

write.table(sub3, file = 'TableS4.txt', quote = F, row.names = F, col.names = F, sep = '\t')
