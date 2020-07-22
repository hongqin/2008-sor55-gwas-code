
rm(list=ls())

files = c("_reformated.3153A_4006873_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv",
"_reformated.Results_SolX_CMI.txt.csv",
"_reformated.Sol125a_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv",
"_reformated.Sol125b_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv",
"_reformated.SolXb_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv",
"_reformated.Sor55_4006777_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv")

wildtype.files = files[c(1,2,5)];   #or control files
mutant.files    = files[c(6,3,4)];  #or experiment files with drugs etc

### read in control data
w1.all = read.delim( wildtype.files[1], header=T, sep="\t", fill=T); 
w2.all = read.delim( wildtype.files[2], header=T, sep="\t", fill=T); 
w3.all = read.delim( wildtype.files[3], header=T, sep="\t", fill=T); 

w1.all = w1.all[ w1.all$geneflag==1, ]
w2.all = w2.all[ w2.all$geneflag==1, ]
w3.all = w3.all[ w3.all$geneflag==1, ]

w1.all = w1.all[ ! is.na(w1.all$geneflag), ]
w2.all = w2.all[ ! is.na(w2.all$geneflag), ]
w3.all = w3.all[ ! is.na(w3.all$geneflag), ]

summary(lm( w1.all$FGMean ~ w1.all$Row + w1.all$Column));
summary(lm( w2.all$FGMean ~ w2.all$Row + w2.all$Column));
summary(lm( w3.all$FGMean ~ w3.all$Row + w3.all$Column));


#rownames(w1.all) = as.character( w1.all$Name ); #some probes are printed >1 times. 

### read in experimental data
m1.all = read.delim( mutant.files[1], header=T, sep="\t", fill=T); 
m2.all = read.delim( mutant.files[2], header=T, sep="\t", fill=T); 
m3.all = read.delim( mutant.files[3], header=T, sep="\t", fill=T); 

m1.all = m1.all[ m1.all$geneflag==1, ];
m2.all = m2.all[ m2.all$geneflag==1, ];
m3.all = m3.all[ m3.all$geneflag==1, ];

m1.all = m1.all[ ! is.na(m1.all$geneflag), ];
m2.all = m2.all[ ! is.na(m2.all$geneflag), ];
m3.all = m3.all[ ! is.na(m3.all$geneflag), ];


summary(lm( m1.all$FGMean ~ m1.all$Row + m1.all$Column));
summary(lm( m2.all$FGMean ~ m2.all$Row + m2.all$Column));
summary(lm( m3.all$FGMean ~ m3.all$Row + m3.all$Column));



### I should normalized them here?
scale = sum(w1.all$FGMean); # w1 is the scale
w2.all$FGMean = w2.all$FGMean * scale / sum(w2.all$FGMean)
w3.all$FGMean = w3.all$FGMean * scale / sum(w3.all$FGMean)
m1.all$FGMean = m1.all$FGMean * scale / sum(m1.all$FGMean)
m2.all$FGMean = m2.all$FGMean * scale / sum(m2.all$FGMean)
m3.all$FGMean = m3.all$FGMean * scale / sum(m3.all$FGMean)

### found out all orfs
orfs = unique( m1.all$orf );
orfs = orfs[ ! is.na(orfs) ]

 # orfs = orfs[1:5];  ############ for debug; 

out = data.frame( matrix(nrow=length(orfs), ncol=7 ) );
names(out) = c( "orf", "chr", "gene", "WtMean","MutMean", "fold.mut.by.wt", "p");

 i = 0;
for( myorf in as.character(orfs) ) {
 i = i + 1;
 wt = c(); mut = c();
 
 v1 = w1.all$FGMean[as.character(w1.all$orf) == myorf ]; v1;
 v2 = w2.all$FGMean[as.character(w2.all$orf) == myorf ]; v2;
 v3 = w3.all$FGMean[as.character(w3.all$orf) == myorf ]; v3;
 wt = c(v1,v2,v3);
 #wt = c(v2,v3);

 u1 = m1.all$FGMean[as.character(m1.all$orf) == myorf ]; u1;
 u2 = m2.all$FGMean[as.character(m2.all$orf) == myorf ]; u2;
 u3 = m3.all$FGMean[as.character(m3.all$orf) == myorf ]; u3;
 mut = c(u1,u2,u3);
 #mut = c( u2, u3);

 t = t.test( wt, mut );
 t2 = t.test( wt, mut, paired=T);

 tmp = w1.all[as.character(w1.all$orf) == myorf , ]
 tmp[1,c("orf","chr","gene")]
 out$chr[i]  = as.character( tmp$chr[1]);
 out$gene[i] = as.character( tmp$gene[1]);
 out$WtMean[i]  = mean( wt );
 out$MutMean[i] = mean( mut );
 out$fold.mut.by.wt[i] = mean(mut) / mean( wt );
 out$orf[i] = myorf;

 out$p[i] = t$p.value; 
 out$p.paired[i] = t2$p.value;
}

 write.table(out, "_out.Sor55.091607.csv", col.names=T, row.names=F, quote=F, sep="\t");

###
 up = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt > 1.5), ]
 write.table(up, "_out.Sor55.up.091707.csv", col.names=T, row.names=F, quote=F, sep="\t");

 down = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt < 0.5), ]
 write.table(down, "_out.Sor55.down.091707.csv", col.names=T, row.names=F, quote=F, sep="\t");


### to do FDR


q("yes");



# average change by chromosome
for( c in 1:7) { #loop over 7 chromosomes
 chr = chrs[c]; 
 tb.sub = tb[tb$Chr == chr,]

 out$avg.change1[c] = mean(tb.sub$f1);
 out$avg.change2[c] = mean(tb.sub$f2);
}


#### 030107v2 beging partition genes by fold using Sol125 only

out2= tb[1,];
out2[1,] = NA;
out2$labels = NA;

for( c in 1:length(chrs) ) {

sub = tb[tb$Chr == chrs[c], ]

q = quantile( sub$f2 , probs=c(0.05,0.1, 0.25, 0.75,0.9, 0.95));
cats = c("low05","low10","low25","avg","up25","up10","up5");
pos  = seq( 1: length(q)-2 );

get.cat <- function( x ) {
 flags   = ifelse( q < x, 1, 0 );
 my.pos = max(pos*flags) + 1;
 if( my.pos > 6 ) { my.pos = 6 }
 ret <- cats[my.pos]
}

 for( i in 1:length(sub[,1]) ) {
   sub$labels[i] = get.cat( sub$f2[i] );
 }

 out2 = rbind(out2, sub);
}

out2 = out2[ ( ! is.na(out2$f2) ), ]

write.table(out2,"030107.partition.by.change.in.chr.csv",row.name=F,col.name=T,sep="\t");

save.image("030107.sorbose.average.RData");

q("no");








