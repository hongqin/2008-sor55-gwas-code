### April 13,version B 
 # additive model 
 #      signal ~ copynum 
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

 ### ANOVA for individual genes
 orfs = as.character(unique( ch5$orf ));

 labels= c("orf","chr");
 out.chr5.FGMean = data.frame(matrix(nrow=length(orfs), ncol=length(labels) ));
 out.chr5.FGMean2 = data.frame(matrix(nrow=length(orfs),ncol=length(labels) ));

 count = 0;

# for( i in 1:10) {
 for( i in 1:length(orfs)) {
    orf = orfs[i];
    sub = allgene[ as.character(allgene$orf)==orf, ]
    out.chr5.FGMean[i,c(1,2)]  = c( orf, 'Chr5');
    out.chr5.FGMean2[i,c(1,2)] = c( orf, 'Chr5');

    m = lm( log2(sub$FGMean) ~ sub$copynum + sub$GC + sub$expt:sub$strain);
    s = summary( m );
    a = anova(m);
    sq = a[,"Sum Sq"];
    out.chr5.FGMean$pf[i] = 1 - pf( s$fstatistic[1], s$fstatistic[2], s$fstatistic[3]);
    out.chr5.FGMean$copynumRole[i] = a[1,2] / sum( sq );
    out.chr5.FGMean$copynumPr[i] = a[1,5];
 
    m = lm( log2(sub$FGMean2) ~ sub$copynum + sub$GC + sub$expt);
    s = summary( m );
    a = anova(m);
    sq = a[,"Sum Sq"];
    out.chr5.FGMean2$pf[i] = 1 - pf( s$fstatistic[1], s$fstatistic[2], s$fstatistic[3]);
    out.chr5.FGMean2$copynumRole[i] = a[1,2] / sum( sq );
    out.chr5.FGMean2$copynumPr[i] = a[1,5];
 
    count = count + 1;
 }


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








q("no");

    # m = lm( log2(sub$FGMean2) ~ sub$copynum );
    # s = summary( m );
    # out.chr5.FGMean2$cpnumRole[i] = s$r.squared[1];


