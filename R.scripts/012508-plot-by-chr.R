### Jan 24, 2008
### plot expression ratio along chromosomes
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

xscale = 3.2E6;  xlim=c(0, xscale); ylim = c(0.1, 4);
cutoffs = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );

pdf( "_012508-chrs.pdf",height=10,width=8);

mat = matrix( seq(1,length(chrs)), nrow=length( chrs), ncol= 2 ); 
layout(mat);

for ( i in 1:length(chrs)) {
  	mychr = chrs[ i ];
	tmp = out[out$chr== mychr,] ;
	exp.chr5 = tmp[ ! is.na(tmp$orf), ]
	
	## the annotations
	tmp = tb[tb$chr == as.character(locs[i]), ];
	chr5 = tmp;
	chr5$orf = as.character( chr5$orf );
	chr5$name = as.character( chr5$name );
	#chr5$len = ifelse( chr5$start > chr5$end, (chr5$start - chr5$end + 1), (chr5$end - chr5$start+1));
	#chr5$flag = ifelse( chr5$len == chr5$cal_len, 1, 0);
	
	## now combine these table
	exp.chr5$start = chr5$start[ match( exp.chr5$orf, chr5$orf ) ];
	exp.chr5$end   = chr5$end  [ match( exp.chr5$orf, chr5$orf ) ]; 
	
	exp.chr5$loc = exp.chr5$start /2  + exp.chr5$end/2 
		
	u = 1; 	#u =  mean( exp.chr5$fold.mut.by.wt); # [1] 0.7891938
	
	par(mar=c(2,2,1,1));
	plot( exp.chr5$fold.mut.by.wt ~ exp.chr5$loc, col="white",xlab='bp',ylab='sor55/wt',log='y',
		xlim=xlim,ylim=ylim); 
	# main= paste( "expression ratio sor55/wt along", mychr, sep=' ') );
	
	for( j in 1:length(exp.chr5$orf) ) {
		yin = exp.chr5$fold.mut.by.wt[j];
		mycol = 'NA';
		if ( yin > cutoffs[2] ) { mycol = 'blue'; 
		} else if ( yin < cutoffs[1] ) { mycol = 'red';
		} else { mycol= 'green'; }
		y = ifelse ( yin > u, yin-u, yin-u );
		MyBox( exp.chr5$loc[j], u, y, xscale/1000, col= mycol);
	}
	
	MyBox( 0, u, 0.01, max(exp.chr5$loc,na.rm=T), col="black");
	points( cens$start[cens$orf==cen.orfs[i]], u+0.005 ,pch=19,cex=0.8, col="black");
	text( -5E4, 2, chrs[i] );
	
}#chr loop

dev.off();


q("no");



cutoffs = quantile( exp.chr5$fold.mut.by.wt, probs=c(0.15,0.85) );
p = c( 0.15, 0.7, 0.15);

for( i in 1:length(exp.chr5$orf) ) {
   fold = exp.chr5$fold.mut.by.wt[i];
   
   type= 99;
   if ( fold > cutoffs[2] ) { type = 'up'; 
   } else if ( fold < cutoffs[1] ) { type = 'down';
   } else { type = 'mean';  }
   
   exp.chr5$type[i] = type;
}

#p = table( exp.chr5$type) / sum(table);
x = c( 30, 70, 0);
pi = exp(loglh.multinomial.sampling( x, p ) ) 
pi;

sorted.locs = sort( exp.chr5$loc );


##sliding window chi-sqaure test?

labels = c('orf','loc', 'type', 'pchisq');
window = 10;
sliding.test = data.frame( matrix(nrow=length(sorted.locs - window + 1), ncol=length(labels) ) );
names( sliding.test) = labels;

sliding.test$loc = sorted.locs;
sliding.test$orf  = exp.chr5$orf [match( sorted.locs, exp.chr5$loc )];
sliding.test$type = exp.chr5$type[match( sorted.locs, exp.chr5$loc )];
sliding.test[1:5,];
exp.chr5[exp.chr5$orf=='orf19.5692',] #check passed.

for( i in 1:length( sliding.test$loc) ) {
   xs = i:(i+window-1);
   types = sliding.test$type[xs];
   downs = length( types[ types=="down"] );
   means = length( types[ types=="mean"] );
   ups = length( types[ types=="up"] );
   obs = c( downs, means, ups);
   expecteds = c(1.5, 7, 1.5);
   k2 = sum( (ups - expecteds)^2 / expecteds );
   sliding.test$pchisq[i] = pchisq( k2, length(obs)-1, lower.tail = F);
}

plot( sliding.test$pchisq ~ sliding.test$loc, log='y');


