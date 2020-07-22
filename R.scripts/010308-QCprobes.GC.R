### use QC probes as normalization 

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

wildtype.files  = files[c(2, 5, 7)];   # control files
mutant.files    = files[c(3, 4, 8)];  # experiment files

### read in control data
w1.all = read.delim( wildtype.files[1], header=T, sep="\t", fill=T); 
w2.all = read.delim( wildtype.files[2], header=T, sep="\t", fill=T); 
w3.all = read.delim( wildtype.files[3], header=T, sep="\t", fill=T); 

w1.all = w1.all[ ! is.na(w1.all$geneflag), ]
w2.all = w2.all[ ! is.na(w2.all$geneflag), ]
w3.all = w3.all[ ! is.na(w3.all$geneflag), ]

pw1.all = w1.all[ w1.all$geneflag==0, ]
pw1.all = pw1.all[ ! is.na(pw1.all$geneflag), ]

library(seqinr)

for( i in 1:length(w1.all[,1]) ) {
 w1.all$GC[i] = GC(s2c(as.character( w1.all$Sequence[i] )) );
 w2.all$GC[i] = GC(s2c(as.character( w2.all$Sequence[i] )) );
 w3.all$GC[i] = GC(s2c(as.character( w3.all$Sequence[i] )) );
}

### read in experimental data
m1.all = read.delim( mutant.files[1], header=T, sep="\t", fill=T); 
m2.all = read.delim( mutant.files[2], header=T, sep="\t", fill=T); 
m3.all = read.delim( mutant.files[3], header=T, sep="\t", fill=T); 

m1.all = m1.all[ ! is.na(m1.all$geneflag), ];
m2.all = m2.all[ ! is.na(m2.all$geneflag), ];
m3.all = m3.all[ ! is.na(m3.all$geneflag), ];

for( i in 1:length(w1.all[,1]) ) {
 m1.all$GC[i] = GC(s2c(as.character( m1.all$Sequence[i] )) );
 m2.all$GC[i] = GC(s2c(as.character( m2.all$Sequence[i] )) );
 m3.all$GC[i] = GC(s2c(as.character( m3.all$Sequence[i] )) );
}


### normalization by QCs

scale = sum( w1.all$FGMean[w1.all$geneflag==1], na.rm=T);

w2.all$FGMean = w2.all$FGMean * scale / sum(w2.all$FGMean[w2.all$geneflag==0], na.rm=T)
w3.all$FGMean = w3.all$FGMean * scale / sum(w3.all$FGMean[w3.all$geneflag==0], na.rm=T)
m1.all$FGMean = m1.all$FGMean * scale / sum(m1.all$FGMean[m1.all$geneflag==0], na.rm=T)
m2.all$FGMean = m2.all$FGMean * scale / sum(m2.all$FGMean[m2.all$geneflag==0], na.rm=T)
m3.all$FGMean = m3.all$FGMean * scale / sum(m3.all$FGMean[m3.all$geneflag==0], na.rm=T)

QCs = unique( pw1.all$orf );

### check GC% effect
summary( lm( w1.all$FGMean ~ w1.all$GC ) );


#stop here. 

save.image("010308.GC.RData");

q("yes");






########################################
###  all genes 

orfs = unique( w1.all$orf );

##  orfs = orfs[1:100]  #debug

labels = c( "orf", "chr", "gene", "geneflag", "WtMean","MutMean", "fold.mut.by.wt", "p");
out = data.frame( matrix(nrow=length( orfs ), ncol=length(labels) ) );
names(out) = labels;

 i = 0;
for( myorf in as.character( orfs ) ) {
 i = i + 1;
 wt = c(); mut = c();
 
 v1 = w1.all$FGMean[as.character(w1.all$orf) == myorf ]; v1;
 v2 = w2.all$FGMean[as.character(w2.all$orf) == myorf ]; v2;
 v3 = w3.all$FGMean[as.character(w3.all$orf) == myorf ]; v3;
 wt = c(v1,v2,v3);

 u1 = m1.all$FGMean[as.character(m1.all$orf) == myorf ]; u1;
 u2 = m2.all$FGMean[as.character(m2.all$orf) == myorf ]; u2;
 u3 = m3.all$FGMean[as.character(m3.all$orf) == myorf ]; u3;
 mut = c(u1,u2,u3);

 t = t.test( wt, mut );
 t2 = t.test( wt, mut, paired=T);

 tmp = w1.all[as.character(w1.all$orf) == myorf , ]
  # tmp[1,c("orf","chr","gene","GC", "geneflag")]
 out$chr[i]  = as.character( tmp$chr[1]);
 out$gene[i] = as.character( tmp$gene[1]);
 out$geneflag[i] = as.character( tmp$geneflag[1]);
 out$GC[i] =  tmp$GC[1];
 out$WtMean[i]  = mean( wt );
 out$MutMean[i] = mean( mut );
 out$fold.mut.by.wt[i] = mean(mut) / mean( wt );
 out$orf[i] = myorf;

 out$p[i] = t$p.value; 
 out$p.paired[i] = t2$p.value;
}

summary(out);
summary(w1.all$FGMean);

out.qc = out[out$geneflag==0, ]

 write.table(out, "_out.Sor55.normalizedbyQC.010308.csv", col.names=T, row.names=F, quote=F, sep="\t");


### level and GC? 


### "non-expressed" genes

 cutoff.wt  =  quantile( out$WtMean,  probs=c(0.05,0.95) );
 cutoff.mut =  quantile( out$MutMean, probs=c(0.05,0.95) );

 bottom = out[ ((out$WtMean < cutoff.wt[1]) & (out$MutMean < cutoff.mut[1])), ]
 tmp = bottom[bottom$geneflag==0,]
 tmp;

save.image( "010308QC.RData" );

### to do FDR


q("yes");


