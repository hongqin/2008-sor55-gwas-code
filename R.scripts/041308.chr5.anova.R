### April 13, anova on chr5 genes, unfinished
### April 12, 2008 ANOVA (regression)

### Feb 20, 08, log2 transformation

### Feb 19, 2008 norm by non-chr5 genes
### background correction by chips, merge results, calculate sor55/wt
### use non-chr5 orfs as normalization, redo with data with re-annotated CSU-ASU genes

rm(list=ls())
 # load( "041208.nonChr5Normalized.anova.RData" );
str(allgene);

ch5 = allgene[ allgene$chr=='Chr5', ]

#####Now use FGMean2
 bk.std = var( all$FGMean[ all$geneflag==0 ] )^0.5;

 sub = ch5[ ch5$FGMean2 > 100, ]  #not good? 

mm4 = lm( log2(sub$FGMean2) ~ sub$copynum + sub$strain + sub$GC + sub$expt );
#summary(mm4);
anova(mm4);

 ## for individual genes
 orfs = as.character(unique( ch5$orf ));

 labels= c("orf","chr","cpnumRole", "cpnumRoleRaw");
 out.chr5 = data.frame(matrix(nrow=length(orfs), 1 ));
 count = 0;
 for( i in 1:length(orfs)) {
    orf = orfs[i];
    sub = allgene[ as.character(allgene$orf)==orf, ]
    out.chr5$orf[i] = orf;
    out.chr5$chr[i] = 'Chr5'; 
    m = lm( log2(sub$FGMean2) ~ sub$copynum );
    s = summary( m );
    out.chr5$cpnumRole[i] = s$r.squared[1];

    m = lm( log2(sub$FGMean) ~ sub$copynum );
    s = summary( m );
    out.chr5$cpnumRoleRaw[i] = s$r.squared[1];

    count = count + 1;
 }





##### Now use ratios

labels = c( "orf","chr","gene","geneflag","FGMean","FGMean2","Column","Row","Sequence","ratio");

rm(all2);
w1.all$ratio = m1.all$FGMean2 / w1.all$FGMean2
#ch5 = cbind( w1.all[ w1.all$chr=='Chr5', labels ], "old1");
all2 = cbind( w1.all[ , labels ], "old1");
names( all2 ) = c(labels, "expt");

w2.all$ratio = m2.all$FGMean2 / w2.all$FGMean2
tmp = cbind( w2.all[, labels ], "old2");
names( tmp ) = c(labels,"expt");
all2 = rbind( all2, tmp );

w3.all$ratio = m3.all$FGMean2 / w3.all$FGMean2
tmp = cbind( w3.all[ , labels ], "new2");
names( tmp ) = c(labels,"expt");
all2 = rbind( all2, tmp );

allgene2 = all2[all2$geneflag==1,]
allgene2$copynum = ifelse( (allgene2$chr=='Chr5') , 'change', 'nonchange');

sub= allgene2[ (allgene2$ratio>0.1)&(allgene2$ratio<5), ]

m = lm( log2(sub$ratio) ~ sub$copynum );
anova(m);




q("yes");


