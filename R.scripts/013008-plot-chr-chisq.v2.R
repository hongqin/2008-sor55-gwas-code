# bug, no cluster box???
### Jan 30, 2008
### sliding window in basepairs
### plot expression ratio along chromosomes overlayed with 
### sliding window chi-square test
### modified from 012508-plot-by-chr.R

rm(list=ls())

 MyBox <- function( x0,y0,h,w,col) {
   xin = c( x0, x0,   x0+w, x0+w );
   yin = c( y0, y0+h, y0+h, y0   );
   polygon(xin,yin, col=col, border=NA);
 }

#  Exactly form of mulinomial sampling 
# vector of x and p
loglh.multinomial.sampling <- function( x, p ) {
  total = sum( x );
  if ( length(x) == length(p) ) { 
        y = p ^ x;
        ret <- lfactorial(total) - sum( lfactorial(x) ) + sum(log(y) );
  } else {
        ret <- NA;
  }
}

################################

#source ("source/011208-Sor55.b.R");
# or 
 load ("011208corrected.RData");

#labels = c( "orf", "chr", "gene", "geneflag", "WtMean","MutMean", "fold.mut.by.wt", "p");
#out = data.frame( matrix(nrow=length( union(orfs, orfs2) ), ncol=length(labels) ) );

## the expression results

cen.orfs = c('CEN1','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7','CENR');
chrs     = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7', 'ChrR');
locs	 = c( '1', '2','3','4','5','6','7', 'R');

tb = read.delim("cal21-annotation.csv");
tb$chr = as.character( tb$chr );
cens = tb[ match(cen.orfs, tb$orf), ]; 

xscale = 3.2E6;  
xlim=c(0, xscale); 
yscale = 4;
ylim = c(0.1, 4);

cutoffs = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );

window.bp = 5E4;  pcutoff = 1E-1;

expecteds.p = c(0.05, 0.9, 0.05);
		
pdf( "_013008-chrs-chisq.pdf",height=10,width=8);

mat = matrix( seq(1,length(chrs)), nrow=length( chrs), ncol= 2 ); 
layout(mat);

#i=4;
for ( i in 1:length(chrs)) {
  	mychr = chrs[ i ];
	tmp = out[out$chr== mychr,] ;
	exp.chrtb = tmp[ ! is.na(tmp$orf), ]
	
	## the annotations
	tmp = tb[tb$chr == as.character(locs[i]), ];
	chrtb = tmp;
	chrtb$orf = as.character( chrtb$orf );
	chrtb$name = as.character( chrtb$name );
	
	## now combine these table
	exp.chrtb$start = chrtb$start[ match( exp.chrtb$orf, chrtb$orf ) ];
	exp.chrtb$end   = chrtb$end  [ match( exp.chrtb$orf, chrtb$orf ) ]; 
	
	exp.chrtb$loc = exp.chrtb$start /2  + exp.chrtb$end/2 
		
	u = 1; 	#u =  mean( exp.chrtb$fold.mut.by.wt); # [1] 0.7891938
	
	par(mar=c(2,2,1,1));
	plot( exp.chrtb$fold.mut.by.wt ~ exp.chrtb$loc, col="white",xlab='bp',ylab='sor55/wt',log='y',
		xlim=xlim,ylim=ylim); 
	# main= paste( "expression ratio sor55/wt along", mychr, sep=' ') );
	
	for( j in 1:length(exp.chrtb$orf) ) {
		yin = exp.chrtb$fold.mut.by.wt[j];
		mycol = 'NA';
		if ( yin > cutoffs[2] ) { mycol = 'blue'; 
		} else if ( yin < cutoffs[1] ) { mycol = 'red';
		} else { mycol= 'green'; }
		
		exp.chrtb$type[j] = mycol;
		
		y = ifelse ( yin > u, yin-u, yin-u );
		MyBox( exp.chrtb$loc[j], u, y, xscale/1000, col= mycol);
	}
	
	MyBox( 0, u, 0.01, max(exp.chrtb$loc,na.rm=T), col="black");
	points( cens$start[cens$orf==cen.orfs[i]], u+0.005 ,pch=19,cex=0.8, col="black");
	text( -5E4, 2, chrs[i] );
	
	#sliding chi-sq test
	sorted.locs = sort( exp.chrtb$loc );
	labels = c('orf','loc', 'type', 'pchisq');
	sliding.test = data.frame( matrix(nrow=length(sorted.locs), ncol=length(labels) ) );
	names( sliding.test) = labels;

	sliding.test$loc  = sorted.locs;
	sliding.test$orf  = exp.chrtb$orf [match( sorted.locs, exp.chrtb$loc )];
	sliding.test$type = exp.chrtb$type[match( sorted.locs, exp.chrtb$loc )];
	sliding.test$clustertype = NA;
	sliding.test[1:5,];
	#exp.chrtb[exp.chrtb$orf=='orf19.5692',] #check passed.

	for( j in 1:(length(sliding.test$loc) - 5) ) {
		start= sliding.test$loc[j];
		subset = sliding.test[ (sliding.test$loc < (start + window.bp ))&(sliding.test$loc >= start), ]
		xs = match( subset$loc, sliding.test$loc);
		window = length( xs );
		#xs = j:(j+window-1);
		types = sliding.test$type[xs];
		downs = length( types[ types=="red"] );
		means = length( types[ types=="green"] );
		ups = length( types[ types=="blue"] );
		obs = c( downs, means, ups);
		obs[is.na(obs)] = 0;
		expecteds = expecteds.p * length(xs);
		k2 = sum( (obs - expecteds)^2 / expecteds );
 		center = floor( mean(xs) ); 
		sliding.test$pchisq[center] = pchisq( k2, length(obs)-1, lower.tail = F);

		yplot= NA;
		if ( sliding.test$pchisq[ center ] < pcutoff ) {
		  if ((obs[1]> expecteds[1]) && (obs[3]< expecteds[3]) ){
		     sliding.test$clustertype[center] = "red"; yplot = 0.15;
		     #rect( sliding.test$loc[j], yplot, sliding.test$loc[j+window-1], yplot+yscale/100,col="red",
		     #     border=NA );
		  } else if ((obs[1] < expecteds[1]) && (obs[3] > expecteds[3]) ) {
		     sliding.test$clustertype[center] = "blue"; yplot = 0.15;
		  } 
		  # MyBox( (sliding.test$loc[center] - xscale/1200),yplot,yscale/100,xscale/600,
		  #		 col=sliding.test$clustertype[center] );
		     rect( sliding.test$loc[j], yplot, sliding.test$loc[max(xs,na.rm=T)], yplot+yscale/10,
		          border=NA,sliding.test$clustertype[center] );
		} else {
			sliding.test$clustertype[center] = "green"; #greens are not shown
		} 
	     rect( sliding.test$loc[j], yplot, sliding.test$loc[max(xs,na.rm=T)], yplot+yscale/10,
	         border=NA,sliding.test$clustertype[center] );
	
	}

	#sliding.test[1:20,]
	write.table(sliding.test, paste("/tmp/__out.Sor55",mychr,"013008.csv", sep="."),
                 col.names=T, row.names=F, quote=F, sep="\t");
	par(new=T);
	plot ( sliding.test$pchisq ~ sliding.test$loc, col='purple', type='l',xlab='',ylab='',log='y', xlim=xlim, axes=F,lty=2); 

}#chr loop

dev.off();


q("no");

