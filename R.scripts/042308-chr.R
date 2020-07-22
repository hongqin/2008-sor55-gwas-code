### 042308 remove hist part. remove x-axis. tick inside

### based 021908-chr-hist.R
### after background correction and non-chr5 nomalization
### histogram of expression ratio along chromosome, then histogram
### modified from 020508 version

rm(list=ls())

 #left ajusted
 MyBox <- function( x0,y0,h,w,col) {
   xin = c( x0, x0,   x0+w, x0+w );
   yin = c( y0, y0+h, y0+h, y0   );
   polygon(xin,yin, col=col, border=NA);
 }

 #centered
 MyBoxC <- function( x0,y0,h,w,col) {
   xin = c( x0-w/2, x0-w/2,   x0+w/2, x0+w/2 );
   yin = c( y0-h/2, y0+h/2,   y0+h/2, y0-h/2   );
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
#source ("source/021908-sor55-norm-merge.b.R");
# or 
# load ("021908.bkgrdCorrected.nonChr5Normalized.RData");
 load ("022008.bkgrdCorrected.nonChr5Normalized.RData");

################# the expression results
cen.orfs = c('CEN1','CENR','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7');
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7');
locs	 = c( '1', 'R','2','3','4','5','6','7');

tb = read.delim("cal21-annotation.csv");
tb$chr = as.character( tb$chr );
cens = tb[ match(cen.orfs, tb$orf), ]; 

xscale = 3.2E6;  
xlim=c(0, xscale); 
yscale = 4;
ylim = c(0.1, 6);
yy = c(0.1, 0.2, 0.5, 1, 2, 4);
xx = pretty( range( tb$end, na.rm=T ) )

cutoffs = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );

#window.bp = 0.5E5;  window = 17;  pcutoff = 1E-6;#too littel blue
#window.bp = 0.5E5;  window = 15;  pcutoff = 1E-6;
#window.bp = 0.5E5;  window = 25;  pcutoff = 1E-4;  #this too loose? 
#window.bp = 0.5E5;  window = 19;  pcutoff = 1E-4;  #this too loose? 
window.bp = 0.5E5;  window = 15;  pcutoff = 1E-4;  #

expecteds = c(0.05, 0.9, 0.05) * window;

 pdf( "_042308-chrs-chisq.pdf",height=9,width=6);
 #jpeg( "_42308-chrs-chisq.jpg",height=2000,width=1500, horizontal=T);
 #jpeg( "_42308-chrs-chisq.jpg",height=900,width=600, horizontal=T);

mat = matrix( seq(1,length(chrs)), nrow=length( chrs), ncol= 2 ); 

layout(mat, heights= c( 1.15, rep(1, nrow(mat)-2), 1.2) );

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

	if ( i == length(chrs) ) {
	  par(mar=c(2,2,0,1));
	  plot( exp.chrtb$fold.mut.by.wt ~ exp.chrtb$loc, col="white",xlab='',ylab='sor55/wt',log='y',
		xlim=xlim,ylim=ylim, axes=F); 
	  axis( 1, at = xx, tcl=0.2, line= 0);
	} else {
	   if ( i == 1) { 
		par(mar=c(0,2,1,1) ); 
	   } else { 
		par(mar=c(0,2,0,1) ); 
	   }
	   plot( exp.chrtb$fold.mut.by.wt ~ exp.chrtb$loc, col="white",xlab='bp',ylab='sor55/wt',log='y',
		xlim=xlim,ylim=ylim, axes=F); 
	  axis( 1, at = xx, labels=F, tcl=0.2);
	} 

	box( );
        axis(2, at= yy, tcl=0.4, line= 0); 
	#axis(2, at=pretty(range(ylim)) ); 

	# main= paste( "expression ratio sor55/wt along", mychr, sep=' ') );
	
	for( j in 1:length(exp.chrtb$orf) ) {
		yin = exp.chrtb$fold.mut.by.wt[j];
		mycol = 'NA';
		if (( yin > cutoffs[2] ) && ( exp.chrtb$WtMean>5)) { mycol = 'blue';   ####remove outliers??!!!!
		} else if ( yin < cutoffs[1] ) { mycol = 'red';
		} else { mycol= 'green'; }
		
		exp.chrtb$type[j] = mycol;
		
		y = ifelse ( yin > u, yin-u, yin-u );
		MyBox( exp.chrtb$loc[j], u, y, xscale/1000, col= mycol);
	}
	
	MyBox( 0, u, 0.01, max(exp.chrtb$loc,na.rm=T), col="black");
	points( cens$start[cens$orf==cen.orfs[i]], u+0.005 ,pch=19,cex=0.8, col="black");
	#text( -5E4, 2, chrs[i] );
	text( -5E3, 4, chrs[i] );
	
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
	exp.chrtb[exp.chrtb$orf=='orf19.5692',] #check passed.

	for( j in 1:(length(sliding.test$loc) - window + 1) ) {
		xs = j:(j+window-1);
		mylocs = sliding.test$loc[xs];  mylocs = mylocs[ ! is.na(mylocs) ];
		
		if (( max(mylocs) - min(mylocs)) < window.bp ) {
		types = sliding.test$type[xs];
		downs = length( types[ types=="red"] );
		means = length( types[ types=="green"] );
		ups = length( types[ types=="blue"] );
		obs = c( downs, means, ups);
		obs[is.na(obs)] = 0;
		#expecteds = c(0.5, 9, 0.5);
		k2 = sum( (obs - expecteds)^2 / expecteds );
 		center = floor( j + window/2 ); 
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
		   MyBoxC( sliding.test$loc[center],yplot,yscale/100,xscale/200,
		  		 col=sliding.test$clustertype[center] );
		  #   rect( sliding.test$loc[j], yplot, sliding.test$loc[j+window-1], yplot+yscale/100,
		  #        border=NA,sliding.test$clustertype[center] );
		} else {
			sliding.test$clustertype[center] = "green"; #greens are not shown
		}
		
		} ## window.bi if loop
 
	}
	#sliding.test[1:20,]

	#par(new=T);
	#plot ( sliding.test$pchisq ~ sliding.test$loc, col='purple', type='l',xlab='',ylab='',log='y', xlim=xlim, axes=F,lty=2); 

}#chr loop
dev.off();


q("no");

