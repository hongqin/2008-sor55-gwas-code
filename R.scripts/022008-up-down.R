### Feb 5, 2008
### parse out up and down genes by chromosomes
### modified from 012908-chr-pie.R

rm(list=ls())

#source ("source/021108-sor55-norm-merge.b.R");

 load ("022008.bkgrdCorrected.nonChr5Normalized.RData");

########################
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr6','Chr4','Chr7','Chr5');

cutoffs = quantile( out$fold.mut.by.wt[out$geneflag==1], probs=c(0.05,0.95) );

out$type = 'middle';

for( i in 1:length( out[,1]) ) {
  if  ( out$fold.mut.by.wt[i] >= cutoffs[2] ) { 
        out$type[i] = 'up05'; 
  } else if ( out$fold.mut.by.wt[i] <= cutoffs[1]) { 
  	out$type[i] = 'bottom05';
  }
}

table( out$type );

write.table(out, "_out.Sor55.up05bottom05.022008.csv", col.names=T, row.names=F, quote=F, sep="\t");

q("no");
