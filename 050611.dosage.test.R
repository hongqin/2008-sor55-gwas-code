# t.test against null hypothesis

# Ch4/7b of the parental strain 3153A is a chimerical product of the reciprocal recombination between what is called Ch4 and Ch7 in the reference strain SC5314.  See Fig. 1B, Fig. 2, and the beginning of the Results. 
# We know exactly, that the brake occurred on Ch4 at approximately 854 kb and at 236 kb of Ch7. Right portion of Ch4 of ~749 kb including the centromere, and left portion of Ch7 of 236 kb were fused creating Ch4/7b.

#If you range the genes of Ch4 and Ch7 according their ORF coordinates from top to bottom, you will see that starting from orf19.5312 down all Ch4 genes will have triploid DNA (aCGH) ratios.  
#These are from 1.3 to 1.7.  Triploid genes are sitting from 854 kb to the right telomere and are fused now with 236 kb of Ch7 starting from left telomere to 236 kb.
 
#We need to answer 2  questions, 1)  What’s an overall transcriptional increase on the trisomic Ch4/7b and 2)  How many genes outside of the monosomic Ch5 and trisomic Ch4/7b significantly changed their expression.

### non4,5,7 normalization
### Sep 29, 09, change negative values to positives; 
### Aug 31-Sep1, 09 update annotation based on cal21-annotation.csv
### Feb 20, 08, log2 transformation for t.test
### Feb 19, 2008 norm by non-chr5 genes
### background correction by chips, merge results, calculate sor55/wt
### use non-chr5 orfs as normalization, redo with data with re-annotated CSU-ASU genes

load( "030411.myannotation.bkgrdCorrected.nonChr457Normalized.RData" );

antb = read.delim("cal21-annotation.csv", colClasses = c("character","character","character","character","integer","integer","integer", "character") );

out$start = antb$start[ match( out$orf, antb$orf ) ]; 

out$ch4.7b = 0;
out$ch4.7b[out$chr=="Chr4" & out$start > out$start[out$orf=='orf19.5312'] ] = 1;
out$ch4.7b[out$chr=="Chr7" & out$start < 236000 ] = 1;

sub4.7b = out[out$ch4.7b==1, ]
summary(sub4.7b$fold.mut.by.wt)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5329  1.0780  1.3620  1.5490  1.7100 46.1600 
#t.test( sub4.7b$fold.mut.by.wt, mu=1.5 ) #trisomy is null 
#p=0.6657 #We cannot reject the null hypothesis that it is trisomy

t.test( log( sub4.7b$fold.mut.by.wt ), mu=log(1.5) ) #trisomy is null 
#t = -4.6274, df = 405, p-value = 4.993e-06

#t.test( sub4.7b$fold.mut.by.wt, mu=1.0 ) #compesation is null
#p=1.6E-6, we can reject the compensation as null
t.test( log(sub4.7b$fold.mut.by.wt), mu=log(1.0) ) #compensation is null
#p< 2.2E-16


sub5 = out[out$chr=='Chr5', ]
summary(sub5$fold.mut.by.wt)
#t.test( sub5$fold.mut.by.wt, mu=1.0) #compensation is null
#p = 5.8E-7
t.test( log( sub5$fold.mut.by.wt), mu=log(1)) #compensation is null
#p< 2.2E-16

#t.test( sub5$fold.mut.by.wt, mu=0.5) #no compensation is null
# p = 2.2E-16
t.test( log(sub5$fold.mut.by.wt), mu=log(0.5) ) #no compensation is null
# p <2.2E-16


sub = out[out$ch4.7b==0 & out$chr != 'Chr5', ]
summary(sub$fold.mut.by.wt)
t.test( sub$fold.mut.by.wt, mu=1.0) #compensation is null

###redo using the empirical distribution as null
t.test( sub4.7b$fold.mut.by.wt, sub$fold.mut.by.wt)
#t = 4.2155, df = 406.76, p-value = 3.073e-05

t.test( sub5$fold.mut.by.wt, sub$fold.mut.by.wt)
#t = -7.2813, df = 531.666, p-value = 1.197e-12

