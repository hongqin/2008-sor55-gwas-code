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

sub = out[out$ch4.7b==0 & out$chr != 'Chr5', ]
summary(sub$fold.mut.by.wt)
summary(sub$fold.mut.by.wt[sub$chr=='Chr4'])
summary(sub$fold.mut.by.wt[sub$chr=='Chr7'])

sub.sig = sub[sub$p.paired<0.01 & (! is.na(sub$orf)), ]
sub.sig = sub.sig[! is.na(sub.sig$orf), ]
summary(sub.sig)
#184 genes


#left  = c(-999,  0.4, 0.6, 0.9, 1.1, 1.9, 2.1);  #032811
#right = c(0.4,   0.6, 0.9, 1.1, 1.9, 2.1, 9999); #032811
left  = c(-999,  0.5, 0.7, 0.9, 1.1, 1.3 );  #032811
right = c(0.5,   0.7, 0.9, 1.1, 1.3, 9999); #032811

report = data.frame( cbind(left, right) );
sub = sub4.7b[sub4.7b$fold.mut.by.wt !=1, ]
for( i in 1:length(left) ) { 
 subB = sub[ (sub$fold.mut.by.wt < right[i])&(sub$fold.mut.by.wt >= left[i] ),  ]; 
 report$freqB[i] = length(subB[,1]) 
}

report$f = report$freqB*100/367
report
#> report
#    left  right freqB         f
#1 -999.0    0.5     0  0.000000
#2    0.5    0.7     7  1.907357
#3    0.7    0.9    21  5.722071
#4    0.9    1.1    42 11.444142
#5    1.1    1.3    77 20.980926
#6    1.3 9999.0   220 59.945504


