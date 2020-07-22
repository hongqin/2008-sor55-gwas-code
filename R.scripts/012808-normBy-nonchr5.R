### Jan 28, 2008
### normalization by non-chr5 chromosomes
### 

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

chrs = c("Chr1","Chr2","Chr3","Chr4","Chr5", "Chr6", "Chr7", "ChrR"); 

### read in control data
w1.all = read.delim( wildtype.files[1], header=T, sep="\t", fill=T); 
w2.all = read.delim( wildtype.files[2], header=T, sep="\t", fill=T); 
w3.all = read.delim( wildtype.files[3], header=T, sep="\t", fill=T); 

w1.all = w1.all[ ! is.na(w1.all$geneflag), ]
w2.all = w2.all[ ! is.na(w2.all$geneflag), ]
w3.all = w3.all[ ! is.na(w3.all$geneflag), ]

### read in experimental data
m1.all = read.delim( mutant.files[1], header=T, sep="\t", fill=T); 
m2.all = read.delim( mutant.files[2], header=T, sep="\t", fill=T); 
m3.all = read.delim( mutant.files[3], header=T, sep="\t", fill=T); 

m1.all = m1.all[ ! is.na(m1.all$geneflag), ];
m2.all = m2.all[ ! is.na(m2.all$geneflag), ];
m3.all = m3.all[ ! is.na(m3.all$geneflag), ];

### normalization by non-chr5 chromosomes
tmp = w1.all[ ! (as.character(w1.all$chr) == "Chr5"), ];
tmp[1:5,c("gene","chr")]
tmp2 = tmp[! is.na(tmp[,1]),]

#scale.old = sum( w1.all$FGMean[ w1.all$geneflag == 1 ], na.rm=T);
scale = sum( w1.all$FGMean[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)], na.rm=T);

w2.all$FGMean = w2.all$FGMean * scale / sum( w2.all$FGMean[(w2.all$chr != 'Chr5') & (w2.all$geneflag==1)], na.rm=T);

w3.all$FGMean = w3.all$FGMean * scale / sum( w3.all$FGMean[(w3.all$chr != 'Chr5') & (w3.all$geneflag==1)], na.rm=T);

m1.all$FGMean = m1.all$FGMean * scale / sum( m1.all$FGMean[(m1.all$chr != 'Chr5') & (m1.all$geneflag==1)], na.rm=T);

m2.all$FGMean = m2.all$FGMean * scale / sum( m2.all$FGMean[(m2.all$chr != 'Chr5') & (m2.all$geneflag==1)], na.rm=T);

m3.all$FGMean = m3.all$FGMean * scale / sum( m3.all$FGMean[(m3.all$chr != 'Chr5') & (m3.all$geneflag==1)], na.rm=T);


########################################
### partition by chromosomes

fields = c(  "chr", "WtMean","MutMean", "fold.mut.by.wt", "p");
out = data.frame( matrix(nrow=length(chrs), ncol= length(fields) ) );
names(out) = fields;

### average change by chromosome
for( c in 1: length(chrs)) { #loop over 7 chromosomes
 chr = chrs[c]; 
 w1a.sub = w1.all$FGMean[ as.character(w1.all$chr) == chr ];
 w2a.sub = w2.all$FGMean[ as.character(w2.all$chr) == chr ];
 w3a.sub = w3.all$FGMean[ as.character(w3.all$chr) == chr ];

 m1a.sub = m1.all$FGMean[ as.character(m1.all$chr) == chr ];
 m2a.sub = m2.all$FGMean[ as.character(m2.all$chr) == chr ];
 m3a.sub = m3.all$FGMean[ as.character(m3.all$chr) == chr ];

 wt  = c( w1a.sub, w2a.sub, w3a.sub);
 mut = c( m1a.sub, m2a.sub, m3a.sub);

 out$chr[c] = chr;
 out$WtMean[c] = mean( wt, na.rm=T );
 out$MutMean[c] = mean( mut, na.rm=T );
 out$fold.mut.by.wt [c] = out$MutMean[c] / out$WtMean[c] 
 
 t = t.test( wt, mut, paired=T);
 out$p[c] = t$p.value;

}

write.table(out, "_out.Sor55.byChr.012808.csv", col.names=T, row.names=F, quote=F, sep="\t");

###

save.image("012808-nonChr5-normalizaiton.RData");

q("no");


