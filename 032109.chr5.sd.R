### 110408, remove 8 genes, remove background level genes
### partition
### 102808 remove some orfs around MTLA locus
###oct11, remove background genes (bottom 5% genes)
###oct10, after fixing chr annotation erros
### 042408 partition chr5 genes

rm(list=ls())
 load( "101008.bkgrdCorrected.nonChr5Normalized.myChr.RData");

# this is the file containing std 
out.sd = read.csv("_out.sor55.nonCh5norm.with.sd.121808.csv");

#for( i in 1:length(out.sd[,1]) ) {
#for( i in 1:10 ) {
# tmp = sqrt( (out.sd$sd.mut[i]/out.sd$MutMean2[i])^2 + (out.sd$sd.wt[i]/out.sd$WtMean2[i])^2 );
 #out.sd$fold.sd[i] = out.sd$fold.mut.by.wt[i] * tmp;
#}

chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr6','Chr4','Chr7','Chr5');

wt.q  =  quantile( out$WtMean2,  probs = c(0.05) );
mut.q =  quantile( out$MutMean2, probs = c(0.05) );

i = 8;
mychr = chrs[ i ];
exp.chr5 = out.sd[out.sd$chr== mychr,] ;
exp.chr5 = exp.chr5[ ! is.na(exp.chr5[,1]), ]

#tb5 = exp.chr5[, c("orf","chr","WtMean","MutMean","WtMean2","MutMean2","fold.mut.by.wt",
#"p","p.paired","sd.wt","sd.mut","sd.mut.vs.wt.ratio")]

tb5 = exp.chr5[,c("orf", "chr", "gene", "WtMean", "MutMean", "WtMean2", "sd.wt", "MutMean2", "sd.mut", "fold.mut.by.wt", "sd.mut.vs.wt.ratio", "p", "p.paired")]

for( j in 4:9 ) {    tb5[,j] = round( tb5[,j], 1) }
for( j in 10:13 ) {    tb5[,j] = round( tb5[,j], 2) }

write.csv( tb5, "_out.chr5.030609.sd.rounded.csv", row.names=F)

quit( "no" );
	
