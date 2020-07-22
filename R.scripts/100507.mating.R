
# Are the mating pathway inhibited in Sor55 mutant.  Sor55 has only 1 copy of Ch5. The Mat
# locus is on Chr5. 

rm(list=ls())

###
BennettTable3 = "Bennett03.table3.csv";
BT3 = read.table( BennettTable3, header=T, sep="\t", fill=T);

myBT3 = BT3[BT3$hqinTag==1, ];   #pick the genes verified by H Qin
orfs = myBT3$ORF19;

### expression data
files = c("_reformated.3153A_4006873_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv",
"_reformated.Results_SolX_CMI.txt.csv",
"_reformated.Sol125a_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv",
"_reformated.Sol125b_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv",
"_reformated.SolXb_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv",
"_reformated.Sor55_4006777_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv")

wildtype.files = files[c(1,2,5)];   ### The orders of the files indicate 
mutant.files    = files[c(6,3,4)];  ### their experimental groups

ceilings = rep(30000,6);

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

 ### todo: correct the bias? 

### I will pick the myBT3$Gene from expression data

myWT = w1.all[1, c("Feature","gene","orf", "Probe","FGMean")]; 
myWT[1,] = t( rep(NA, length(myWT[1,])) );
myMT = myWT; 

i = 0;
 # tmp = as.character(orfs);  ### for debug
 # myorf = tmp[i+1];            ### for debug

#for( myorf in as.character(orfs[1:3]) ) {
for( myorf in as.character(orfs) ) {
 i = i + 1;

 v1 = w1.all[as.character(w1.all$orf) == myorf,  ]; v1[, c("ID","orf", "Probe")]
 v2 = w2.all[as.character(w2.all$orf) == myorf,  ]; v2[, c("ID","orf", "Probe")]
 v3 = w3.all[as.character(w3.all$orf) == myorf,  ]; v3[, c("ID","orf", "Probe")]
 mytmp = rbind(v1[, c("Feature","gene","orf", "Probe","FGMean")], 
	v2[, c("Feature","gene","orf", "Probe","FGMean")],
	v3[, c("Feature","gene","orf", "Probe","FGMean")]);
 myWT = rbind( myWT, mytmp);

 u1 = m1.all[as.character(m1.all$orf) == myorf,  ]; u1[, c("ID","orf", "Probe")]
 u2 = m2.all[as.character(m2.all$orf) == myorf,  ]; u2[, c("ID","orf", "Probe")]
 u3 = m3.all[as.character(m3.all$orf) == myorf,  ]; u3[, c("ID","orf", "Probe")]
 mytmp2 = rbind( u1[, c("Feature","gene","orf", "Probe","FGMean")], 
	u2[, c("Feature","gene","orf", "Probe","FGMean")],
	u3[, c("Feature","gene","orf", "Probe","FGMean")]);
 myMT = rbind( myMT, mytmp2);
}

 myWT = myWT[ ! is.na(myWT$Probe), ]  #Remove the first NA row
 myMT = myMT[ ! is.na(myMT$Probe), ]  #Remove the first NA row

### check the list; 
 ii = ceiling( runif(5)*100);
 myWT[ ii, c("orf","Probe")]
 myMT[ ii, c("orf","Probe")]
 as.logical( myWT$Probe==myMT$Probe );

### pairwise t-test; 
 t.test( myWT$FGMean, myMT$FGMean, paired=T);  # p-value = 0.001038 YES. 
 mean(myMT$FGMean) / mean(myWT$FGMean)  #80% increase 

quit("yes"); 

###############################
###
 up = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt > 1.5), ]
 write.table(up, "_out.Sor55.up.091707.csv", col.names=T, row.names=F, quote=F, sep="\t");

 down = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt < 0.5), ]
 write.table(down, "_out.Sor55.down.091707.csv", col.names=T, row.names=F, quote=F, sep="\t");

### to do FDR

#### 030107v2 beging partition genes by fold using Sol125 only

out2= tb[1,];
out2[1,] = NA;
out2$labels = NA;

for( c in 1:length(chrs) ) {
sub = tb[tb$Chr == chrs[c], ]; 

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








