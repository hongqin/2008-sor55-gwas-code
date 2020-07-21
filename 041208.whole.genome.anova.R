### April 12, 2008 ANOVA (regression)

### Feb 20, 08, log2 transformation?? aborted afterwards

### Feb 19, 2008 norm by non-chr5 genes
### background correction by chips, merge results, calculate sor55/wt
### use non-chr5 orfs as normalization, redo with data with re-annotated CSU-ASU genes

rm(list=ls())
files = c(
"_reformated.3153A_4006873_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv",  #1 wt uneven
"_reformated.Results_SolX_CMI.txt.csv",		#2 wt
"_reformated.Sol125a_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #3 mut
"_reformated.Sol125b_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #4 mut
"_reformated.SolXb_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #5  wt
"_reformated.Sor55_4006777_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv", #6  mut uneven
"_reformated.3153A-121407.csv.csv", # 7 wt
"_reformated.Sor55-121407.csv.csv"  #8  mut
);

wildtype.files  = files[c(2, 5, 7, 1)];   # control files
mutant.files    = files[c(3, 4, 8, 6)];  # experiment files

### read in control data
w1.all = read.delim( wildtype.files[1], header=T, sep="\t", fill=T); 
w2.all = read.delim( wildtype.files[2], header=T, sep="\t", fill=T); 
w3.all = read.delim( wildtype.files[3], header=T, sep="\t", fill=T); 
w4.all = read.delim( wildtype.files[4], header=T, sep="\t", fill=T); 

w1.all = w1.all[ ! is.na(w1.all$geneflag), ]
w2.all = w2.all[ ! is.na(w2.all$geneflag), ]
w3.all = w3.all[ ! is.na(w3.all$geneflag), ]
w4.all = w4.all[ ! is.na(w4.all$geneflag), ]

#pw1.all = w1.all[ w1.all$geneflag==0, ]
#pw2.all = w2.all[ w2.all$geneflag==0, ]
#pw3.all = w3.all[ w3.all$geneflag==0, ]
#pw1.all = pw1.all[ ! is.na(pw1.all$geneflag), ]
#pw2.all = pw2.all[ ! is.na(pw2.all$geneflag), ]
#pw3.all = pw3.all[ ! is.na(pw3.all$geneflag), ]

### read in experimental data
m1.all = read.delim( mutant.files[1], header=T, sep="\t", fill=T); 
m2.all = read.delim( mutant.files[2], header=T, sep="\t", fill=T); 
m3.all = read.delim( mutant.files[3], header=T, sep="\t", fill=T); 
m4.all = read.delim( mutant.files[4], header=T, sep="\t", fill=T); 

m1.all = m1.all[ ! is.na(m1.all$geneflag), ];
m2.all = m2.all[ ! is.na(m2.all$geneflag), ];
m3.all = m3.all[ ! is.na(m3.all$geneflag), ];
m4.all = m3.all[ ! is.na(m4.all$geneflag), ];

#pm1.all = m1.all[ m1.all$geneflag==0, ];
#pm2.all = m2.all[ m2.all$geneflag==0, ];
#pm3.all = m3.all[ m3.all$geneflag==0, ];
#pm1.all = pm1.all[ ! is.na(pm1.all$geneflag), ];
#pm2.all = pm2.all[ ! is.na(pm2.all$geneflag), ];
#pm3.all = pm3.all[ ! is.na(pm3.all$geneflag), ];

### background corection

w1.background   = median( w1.all$FGMean[w1.all$geneflag==0])
w2.background   = median( w2.all$FGMean[w2.all$geneflag==0])
w3.background   = median( w3.all$FGMean[w3.all$geneflag==0])
w4.background   = median( w4.all$FGMean[w4.all$geneflag==0])
c(w1.background,w2.background,w3.background ,w4.background )

m1.background   = median( m1.all$FGMean[m1.all$geneflag==0])
m2.background   = median( m2.all$FGMean[m2.all$geneflag==0])
m3.background   = median( m3.all$FGMean[m3.all$geneflag==0])
m4.background   = median( m4.all$FGMean[m4.all$geneflag==0])
c(m1.background,m2.background,m3.background ,m4.background )

w1.all$FGMean2 = w1.all$FGMean - w1.background;
w2.all$FGMean2 = w2.all$FGMean - w2.background;
w3.all$FGMean2 = w3.all$FGMean - w3.background;
w4.all$FGMean2 = w4.all$FGMean - w4.background;

m1.all$FGMean2 = m1.all$FGMean - m1.background;
m2.all$FGMean2 = m2.all$FGMean - m2.background;
m3.all$FGMean2 = m3.all$FGMean - m3.background;
m4.all$FGMean2 = m4.all$FGMean - m4.background;

### normalization

#scale = sum( w1.all$FGMean2[w1.all$geneflag==1], na.rm=T);
scale = sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)], na.rm=T);

w2.all$FGMean2 = w2.all$FGMean * scale / sum( w2.all$FGMean2[(w2.all$chr != 'Chr5') &   (w2.all$geneflag==1)],na.rm=T);

w3.all$FGMean2 = w3.all$FGMean * scale / sum( w3.all$FGMean2[(w3.all$chr != 'Chr5') & (w3.all$geneflag==1)],na.rm=T);

w4.all$FGMean3 = w4.all$FGMean * scale / sum( w4.all$FGMean2[(w4.all$chr != 'Chr5') & (w4.all$geneflag==1)],na.rm=T);

m1.all$FGMean2 = m1.all$FGMean * scale / sum( m1.all$FGMean2[(m1.all$chr != 'Chr5') & (m1.all$geneflag==1)],na.rm=T);

m2.all$FGMean2 = m2.all$FGMean * scale / sum( m2.all$FGMean2[(m2.all$chr != 'Chr5') & (m2.all$geneflag==1)],na.rm=T);

m3.all$FGMean2 = m3.all$FGMean * scale / sum( m3.all$FGMean2[(m3.all$chr != 'Chr5') & (m3.all$geneflag==1)],na.rm=T);

m4.all$FGMean2 = m4.all$FGMean * scale / sum( m4.all$FGMean2[(m4.all$chr != 'Chr5') & (m4.all$geneflag==1)],na.rm=T);

### do some checks
w3.all[w3.all$orf=="orf19.2896",c("FGMean","FGMean2")] #SOU1 has 12 probes


########################################
### found out all genes 

orfs  = as.character( unique( w1.all$orf ));
orfs2 = as.character( unique( w3.all$orf ));
orfs.new = setdiff( orfs2, orfs);
tmp = union( orfs, orfs2);


#wildtype.files  = files[c(2, 5, 7, 1)];   # control files
#mutant.files    = files[c(3, 4, 8, 6)];  # experiment files

labels = c( "orf","chr","gene","geneflag","FGMean","FGMean2","Column","Row","Sequence");

rm( all );
all = cbind( w1.all[, labels ], "WT", "old1");
names( all ) = c(labels, "strain", "expt");

tmp = cbind( w2.all[, labels ], "WT", "old2");
names( tmp ) = c(labels, "strain", "expt");
all = rbind( all, tmp );

tmp = cbind( w3.all[, labels ], "WT", "new2");
names( tmp ) = c(labels, "strain", "expt");
all = rbind( all, tmp );

tmp = cbind( m1.all[, labels ], "Sor55", "old1");
names( tmp ) = c(labels, "strain", "expt");
all = rbind( all, tmp );

tmp = cbind( m2.all[, labels ], "Sor55", "old2");
names( tmp ) = c(labels, "strain", "expt");
all = rbind( all, tmp );

tmp = cbind( m3.all[, labels ], "Sor55", "new2");
names( tmp ) = c(labels, "strain", "expt");
all = rbind( all, tmp );

all$copynum = ifelse( ((all$chr=='Chr5') & (all$strain=='Sor55')), 1, 2);

library(seqinr);
#for( i in 1:100 ){
for( i in 1:length(all[,1]) ){ #this takes 150 min on laptop
  if ( is.na(all$Sequence[i]) ) {
    all$GC[i] = NA;
    all$len[i] = NA;
  } else {
    all$GC[i] = GC( s2c( as.character(all$Sequence[i]) ) );
    all$len[i] = length( s2c( as.character(all$Sequence[i]) ) ); 
  }
}

#i am here
allgene = all[ all$geneflag == 1, ]


save.image( "041208.nonChr5Normalized.anova.RData" );

###### use FGMean
#mx = lm( log2(allgene$FGMean) ~ allgene$chr );

m3a = lm( log2(allgene$FGMean) ~ allgene$copynum + allgene$strain + allgene$GC );
summary(m3a);
anova(m3a); #why anova and summary give different p on copynum?
 #anova use F, summary use t

m4 = lm( log2(allgene$FGMean) ~ allgene$copynum + allgene$strain + allgene$Column + allgene$Row );
summary(m4);  ##Column is not significant
anova(m4);

m4b = lm( log2(allgene$FGMean) ~ allgene$copynum + allgene$strain + allgene$GC + allgene$Row);
summary(m4b);
anova(m4b);

#####Now use FGMean2
 bk.std = var( all$FGMean[ all$geneflag==0 ] )^0.5;

 #sub = allgene[ allgene$FGMean2 > (2*bk.std), ]  #not good? 
 #str(sub); # 35537 obs. of  14 variables:

 #sub = allgene[ allgene$FGMean2 > (1*bk.std), ]  #not good? 
 #str(sub); #53585 obs. of  14 variables:

 sub = allgene[ allgene$FGMean2 > 100, ]  #not good? 

#mm4 = lm( log2(sub$FGMean2) ~ sub$copynum + sub$strain + sub$GC + sub$Row ); #Column is not significant 
mm4 = lm( log2(sub$FGMean2) ~ sub$copynum + sub$strain + sub$GC + sub$expt ); #Column is not significant 
#summary(mm4);
anova(mm4);

mm3 = lm( log2(sub$FGMean2) ~ sub$copynum + sub$strain + sub$GC ); 
summary(mm3);
anova(mm3);

mm1 = lm( log2(sub$FGMean2) ~ sub$strain); 
summary(mm1);
anova(mm1);

mm1 = lm( log2(sub$FGMean2) ~ sub$copynum); 
summary(mm1);
anova(mm1);

##### Now use ratios



q("yes");


