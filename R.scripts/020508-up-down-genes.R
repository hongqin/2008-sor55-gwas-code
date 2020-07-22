### Feb 5, 2008
### parse out up and down genes by chromosomes
### modified from 012908-chr-pie.R

rm(list=ls())

#source ("source/011208-Sor55.b.R");
# or 
 load ("011208corrected.RData");

#################### background correction here
out.backup = out;

wt.background   = median(out$WtMean[out$geneflag==0])
mut.background = median(out$MutMean[out$geneflag==0])

out$WtMean.adj  = out$WtMean  - wt.background;
out$MutMean.adj = out$MutMean - mut.background;

 out$WtMean.adj[out$WtMean.adj<0] = 1;
out$MutMean.adj[out$MutMean.adj<0] = 1;

out$fold2.mut.by.wt = out$MutMean.adj / out$WtMean.adj; 


plot(out$fold2.mut.by.wt ~ out$fold.mut.by.wt, xlab='raw ratio',ylab='adj ratio', main ='all data');

#now remove outliers
plot(out$fold2.mut.by.wt[out$WtMean.adj>5] ~ out$fold.mut.by.wt[out$WtMean.adj>5], 
xlab='raw ratio',ylab='adj ratio', main ='low valued removed');

#now only outliers
plot(out$fold2.mut.by.wt[out$WtMean.adj<5] ~ out$fold.mut.by.wt[out$WtMean.adj<5], 
xlab='raw ratio',ylab='adj ratio', main ='low valued removed');


########################
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr6','Chr4','Chr7','Chr5');

cutoffs = quantile( out$fold2.mut.by.wt, probs=c(0.05,0.95) );

out$type = 'middle';

for( i in 1:length( out[,1]) ) {
  if ( ( out$fold2.mut.by.wt[i] > cutoffs[2] )&&(out$WtMean.adj[i]>5) ){ 
        out$type[i] = 'up05'; 
  } else if (( out$fold2.mut.by.wt[i] < cutoffs[1]) && ( out$WtMean.adj[i]>5)) { 
  	out$type[i] = 'bottom05';
  }
}

write.table(out, "_out.Sor55.up05bottom05.020508.csv", col.names=T, row.names=F, quote=F, sep="\t");

q("no");