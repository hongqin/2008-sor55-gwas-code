###oct11, remove background genes (bottom 5% genes)
###oct10, after fixing chr annotation erros
### 042408 partition chr5 genes

rm(list=ls())
 load( "101008.bkgrdCorrected.nonChr5Normalized.myChr.RData");

chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr6','Chr4','Chr7','Chr5');

wt.q  =  quantile( out$WtMean2,  probs = c(0.05) );
mut.q =  quantile( out$MutMean2, probs = c(0.05) );

#cutoffs.fold = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );
#cutoffs.mut  = quantile( out$MutMean2, probs=c(0.25,0.75) );

i = 8;
mychr = chrs[ i ];
tmp = out[out$chr== mychr,] ;
exp.chr5 = tmp[ ! is.na(tmp$orf), ] #507 genes

tmp = exp.chr5[ exp.chr5$WtMean2 < wt.q[1] | exp.chr5$MutMean2 < mut.q[1] , ]

exp.chr5B = exp.chr5[ exp.chr5$WtMean2 > wt.q[1] & exp.chr5$MutMean2 > mut.q[1] , ]
#459 genes

exp.chr5 = exp.chr5B; 

left = seq( 0.0, 2, by=0.1); #left boundries
right = left + 0.1;  #right boundries
left[1] = -999999;
right[length(right)] = 9999;
out = data.frame( cbind(left, right) );

for( i in 1:length(left) ) {
 sub = exp.chr5[ (exp.chr5$fold.mut.by.wt < right[i])&(exp.chr5$fold.mut.by.wt >= left[i] ),  ]; 
 out$freq[i] = length(sub[,1])

}



quit( "no" );
	


