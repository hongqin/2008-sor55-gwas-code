### Jan 29, 2008
### plot expression ratio along chromosomes
rm(list=ls())

#source ("source/011208-Sor55.b.R");
# or 
 load ("011208corrected.RData");

# chrs     = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7', 'ChrR');
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr6','Chr4','Chr7','Chr5');

cutoffs = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );

pdf( "_012908-chr-pie.pdf",height=2,width=4);

mat = matrix( seq(1,10), nrow=2, ncol= 5, byrow=T );
#heights = c(1.2, 1,1,1, 1.2, 1,1,1); 
layout(mat);
# layout.show(8);

for ( i in 1:length(chrs)) {
  	mychr = chrs[ i ];
	tmp = out[out$chr== mychr,] ;
	exp.chr5 = tmp[ ! is.na(tmp$orf), ]
	
	par(mar=c(0.2,0.2,1,0.2));

	frac = NA;
	for ( j in 1:length(exp.chr5$orf) ) {
		yin = exp.chr5$fold.mut.by.wt[j];
		mycol = 'NA';
		if ( yin > cutoffs[2] ) { mycol = 'blue'; 
		} else if ( yin < cutoffs[1] ) { mycol = 'red';
		} else { mycol= 'green'; }
		frac = c( frac, mycol);
		frac = frac[ ! is.na(frac) ];
	}
	dist = table( frac );
	pie( dist, main= mychr, labels=NA, col=c('blue','green','red'),angle=90 );
	
}#chr loop

dev.off();


q("no");


