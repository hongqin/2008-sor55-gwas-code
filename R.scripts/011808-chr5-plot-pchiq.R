### Jan 18, 2008
### plot expression ratio along chr5 
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
tmp = out[out$chr=='Chr5',] ;
exp.chr5 = tmp[ ! is.na(tmp$orf), ]

## the annotations
tb = read.delim("cal21-annotation.csv");

tmp = tb[tb$chr == '5', ];
chr5 = tmp;
chr5$orf = as.character( chr5$orf );
chr5$name = as.character( chr5$name );

chr5$len = ifelse( chr5$start > chr5$end, (chr5$start - chr5$end + 1), (chr5$end - chr5$start+1));
chr5$flag = ifelse( chr5$len == chr5$cal_len, 1, 0);

cen.orfs = c('CEN1','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7','CENR');
cens = tb[ match(cen.orfs, tb$orf), ]; 

## now combine these table
str( exp.chr5$orf );
str( chr5$orf );

match( exp.chr5$orf, chr5$orf);

#exp.chr5$start = chr5$start[ match( chr5$orf , exp.chr5$orf), ];wrong
#exp.chr5$end   = chr5$end  [ match( chr5$orf , exp.chr5$orf), ];wrong
exp.chr5$start = chr5$start[ match( exp.chr5$orf, chr5$orf ) ];
exp.chr5$end   = chr5$end  [ match( exp.chr5$orf, chr5$orf ) ]; 

exp.chr5$loc = exp.chr5$start /2  + exp.chr5$end/2 

# visual check ... OK
 chr5[ chr5$orf=='orf19.4257',]
 exp.chr5[ exp.chr5$orf=='orf19.4257',]

 chr5[ chr5$orf=='orf19.4283',]
 exp.chr5[ exp.chr5$orf=='orf19.4283',]

############################plot  chr5

u =  mean( exp.chr5$fold.mut.by.wt); # [1] 0.7891938
xscale = max(exp.chr5$loc,na.rm=T);

pdf("012308-chr5.pdf",height=3,width=5);
plot( exp.chr5$fold.mut.by.wt ~ exp.chr5$loc, col="white",xlab='bp',ylab='sor55/wt', log='y' );
MyBox( 0, u, 0.01, max(exp.chr5$loc,na.rm=T), col="black");
points( cens$start[cens$orf=='CEN5'], u+0.005 ,pch=19,cex=0.8, col="black");

for( i in 1:length(exp.chr5$orf) ) {
 yin = exp.chr5$fold.mut.by.wt[i];
 mycol = ifelse( yin>u, 'blue','red' );
 #y = ifelse ( yin > u, yin-u+0.005, yin-u-0.005 );
 y = ifelse ( yin > u, yin-u, yin-u );
 #y = yin - u;
 MyBox( exp.chr5$loc[i], u, y, xscale/800, col= mycol);
}
dev.off();




cutoffs = quantile( exp.chr5$fold.mut.by.wt, probs=c(0.15,0.85) );

pdf("012308-chr5-v2.pdf",height=3,width=5);
plot( exp.chr5$fold.mut.by.wt ~ exp.chr5$loc, col="white",xlab='bp',ylab='sor55/wt', log='y' );
MyBox( 0, 2, 0.01, max(exp.chr5$loc,na.rm=T), col="black");
points( cens$start[cens$orf=='CEN5'], 2+0.005 ,pch=19,cex=0.8, col="black");

for( i in 1:length(exp.chr5$orf) ) {
 yin = exp.chr5$fold.mut.by.wt[i];
 mycol = 'NA';
   if ( yin > cutoffs[2] ) { mycol = 'blue'; 
   } else if ( yin < cutoffs[1] ) { mycol = 'red';
   } else { mycol= 'gray'; }
 #y = ifelse ( yin > u, yin-u+0.005, yin-u-0.005 );
 #y = ifelse ( yin > u, yin-u, yin-u );
 y = yin;
 MyBox( exp.chr5$loc[i], 0.001, y, xscale/800, col= mycol);
}
dev.off();






###end of plot of chr5


hist( exp.chr5$fold.mut.by.wt, brk=20 );

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

## to do here. 

###





#
#
q("no");
#
#


###

 up = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt > 1.5), ]
 down = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt < 0.5), ]

 out.qc = out[out$geneflag==0, ]
 cutoff.wt  =  quantile( out$WtMean,  probs=c(0.05,0.95) );
 cutoff.mut =  quantile( out$MutMean, probs=c(0.05,0.95) );

 bottom = out[ ((out$WtMean < cutoff.wt[1]) & (out$MutMean < cutoff.mut[1])), ]
 tmp = bottom[bottom$geneflag==0,]
 tmp; #there are 25 control here, 


