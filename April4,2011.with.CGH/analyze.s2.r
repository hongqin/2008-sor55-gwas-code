s2 =read.csv( "s2.csv")
#names(s2) = c("orf","gene","chr","WtMean","MutMean", "fold.mut.by.wt", "sd", "p", "p.paired", "CGH", "sd2", "start","end", "strand")

x = table( s2$fold.mut.by.wt)
max(x) #726 genes
max(x)/6001 #12.1%

s2$ch4.7b = 0;
s2[s2$orf=='orf19.5312',]
s2$ch4.7b[s2$A21_CHROMOSOME=="Ca21chr4" & s2$A21_START > s2$A21_START[s2$orf=='orf19.5312'] ] = 1;
s2$ch4.7b[s2$chr=="Ca21chr7" & s2$A21_START < 236000 ] = 1;

sub4.7b = s2[s2$ch4.7b==1, ]
summary(sub4.7b$fold.mut.by.wt)
#> summary(sub4.7b$fold.mut.by.wt)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5329  1.0520  1.3330  1.4170  1.6750  4.3590 

sub = s2[s2$ch4.7b==0 & s2$A21_CHROMOSOME != 'Ca21chr5', ] #5208 genes

summary(sub$fold.mut.by.wt)
summary(sub$fold.mut.by.wt[sub$A21_CHROMOSOME =='Ca21chr4'])
summary(sub$fold.mut.by.wt[sub$A21_CHROMOSOME =='Ca21chr7'])

sub.sig = sub[sub$p.paired<0.01 & (! is.na(sub$orf)), ]
sub.sig = sub.sig[! is.na(sub.sig$orf), ]
summary(sub.sig)
#203 genes 




