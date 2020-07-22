q("no");

### April 13,version C, does not work yet.

 # multiplicative model    
 #      sor55Signal ~ WTSignal * (1/copynum);

### anova on chr5, gene by gene method

### April 12, 2008 ANOVA (regression)

### Feb 20, 08, log2 transformation

### Feb 19, 2008 norm by non-chr5 genes
### background correction by chips, merge results, calculate sor55/wt
### use non-chr5 orfs as normalization, redo with data with re-annotated CSU-ASU genes

rm(list=ls())
  load( "041208.nonChr5Normalized.anova.RData" );
str(allgene);

ch5 = allgene[ allgene$chr=='Chr5', ]
#orfs = as.character(unique( ch5$orf[ch5$expt!='new2'] ));

wt = ch5[ ch5$strain=='WT',  ]
sor55 = ch5[ ch5$strain=='Sor55',  ]

ch5.2 = cbind( wt, sor55 );
  #visual check
  ch5.2[1:10,c(1,7:8, 21:22)]

label1 = names(wt);
label2 = paste( label1, 'Sor55', sep='');
names(ch5.2) = c( label1, label2);

ch5.2$FGMean2[ch5.2$FGMean2 < 100] = NA;
ch5.2$FGMean2Sor55[ch5.2$FGMean2Sor55 < 100] = NA;

#m =  lm( log2(ch5.2$FGMean2) ~ log2(ch5.2$FGMean2Sor55) );
#m =  lm( log2(ch5.2$FGMean2) ~ log2(ch5.2$FGMean2Sor55) + ch5.2$GC);
#m =  lm( log2(ch5.2$FGMean2) ~ log2(ch5.2$FGMean2Sor55) + ch5.2$GC + ch5.2$expt );
#summary( lm( log2(ch5.2$FGMean2) ~ log2(ch5.2$FGMean2Sor55) + ch5.2$GC ));
#anova( lm( log2(ch5.2$FGMean2) ~ log2(ch5.2$FGMean2Sor55) + ch5.2$GC ));

m =  lm( log2(ch5.2$FGMean2Sor55) ~ log2(ch5.2$FGMean2) + ch5.2$GC);
summary(m);
anova(m);

plot( log2(ch5.2$FGMean2) ~ log2(ch5.2$FGMean2Sor55) );

plot( log2(ch5.2$FGMean2) ~ ch5.2$GC );

### simple model on scale by MLE
## 4pm, I should use MLE to calculate the scale.
##  sor55singal ~ WTsignal * scale + GC * b + error
## simple model argues for scale = 0.5 because parental strain has 2 copies, 
## while sor55 has 1 copy of chr5

ch5.3 =  data.frame( ch5.2$"orf" );
ch5.3$sor55= log2(ch5.2$FGMeanSor55);
ch5.3$wt =   log2(ch5.2$FGMean2) ;
ch5.3$GC =   ch5.2$GC ;

ch5.3$sor55[ ch5.3$sor55<7] = NA;
ch5.3$wt[ ch5.3$wt<7] = NA;


m = lm( ch5.3$sor55 ~ ch5.3$wt );
plot( ch5.3$sor55 ~ ch5.3$wt );
abline( m, col="red");

m = lm( ch5.3$sor55 ~ ch5.3$wt + ch5.3$GC);
anova(m);
summary(m);



## to do:MLE
## sor55 = wt * (scale <- N(u,sigma) )
##  ie  log(sor55/wt) ~ N(u,sigma);

w1.all$ratio = m1.all$FGMean2 / w1.all$FGMean2;
w2.all$ratio = m2.all$FGMean2 / w2.all$FGMean2;
w3.all$ratio = m3.all$FGMean2 / w3.all$FGMean2;

tmp1 = w1.all[ w1.all$chr=='Chr5', ]
tmp2 = w2.all[ w1.all$chr=='Chr5', ]
tmp3 = w3.all[ w1.all$chr=='Chr5', ]

ch5 = rbind( tmp1[,c("orf","chr","ratio")], tmp2[,c("orf","chr","ratio")], tmp3[,c("orf","chr","ratio")]);

cutoffs = quantile( ch5$ratio, probs = c(0.02,  0.98) );
ch5$ratio[ch5$ratio<cutoffs[1]] = NA; # I need to use quantile as cutoff
ch5$ratio[ch5$ratio>cutoffs[2]] = NA;  # need quantile cutoff

ch5 = ch5[ ! is.na(ch5$orf), ]
ch5 = ch5[ ! is.na(ch5$ratio), ]


hist(ch5$ratio);
summary(ch5$ratio);
# without cutoffs
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-232.3000    0.4518    0.6410    0.8919    0.9952  325.9000 

#> summary(ch5$ratio); Using 2% cutoffs
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1165  0.5010  0.6743  0.9028  1.0160  6.2790 

u = seq( -5, 5, by=0.1);
sigma = seq(0.1, 2, 0.1);
llikm = matrix(nrow=length(u), ncol=length(sigma) ) ;
ztmp = matrix(nrow=length(u), ncol=length(sigma) ) ;

for(i in 1:length(u) ) {
 for( j in 1:length(sigma) ) {
    z = log10( dnorm( log2(ch5$ratio), u[i], sigma[j] ) );
    z[z== -Inf] = NA;
    llikm[i,j] = sum(  z, na.rm=T );
    ztmp[i,j] = i;
 }
}

image(u,sigma,ztmp, col=terrain.colors(10));


image(u,sigma,llikm, col=terrain.colors(10));





q("no");

 
##simulate the role of copynum
N=500;
role = 1;
x2 = rnorm(N) * 10;
x1 = x2 * role * 1/2 + rnorm(N);

x = c(x1,x2);
copynum = c(rep(1,N), rep(2,N));

y = x + rnorm(N);

summary( lm( y ~ copynum ) ); ##Not good


#########try multiplicative model

 to do in version c.







    # m = lm( log2(sub$FGMean2) ~ sub$copynum );
    # s = summary( m );
    # out.chr5.FGMean2$cpnumRole[i] = s$r.squared[1];


