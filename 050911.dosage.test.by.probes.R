### non4,5,7 normalization
### Sep 29, 09, change negative values to positives; 
### Aug 31-Sep1, 09 update annotation based on cal21-annotation.csv
### Feb 20, 08, log2 transformation for t.test
### Feb 19, 2008 norm by non-chr5 genes
### background correction by chips, merge results, calculate sor55/wt
### use non-chr5 orfs as normalization, redo with data with re-annotated CSU-ASU genes

rm(list=ls())

load( "030411.myannotation.bkgrdCorrected.nonChr457Normalized.RData" );

### found out all genes 
#orfs  = as.character( unique( w1.all$orf[ (! is.na(w1.all$geneflag)) & ( w1.all$geneflag == 1) ] ));
#orfs = orfs[ -grep("wrong", orfs) ] #090109 remove the wrong orf
#orfs is in loaded image

#labels = c( "orf", "chr", "gene", "geneflag", "WtMean","MutMean","WtMean2","MutMean2", "sd.wt","sd.mut","fold.mut.by.wt", "sd.ratio","p", "p.paired");
#out = data.frame( matrix(nrow=length( orfs ), ncol=length(labels) ) );
#out = data.frame( matrix(nrow=length( union(orfs, orfs2) ), ncol=length(labels) ) );
#names(out) = labels;
#out is in the loaded image

yftb = read.csv("ch47b.csv"); #tb from Yang Feng

orfs47b = yftb$ORF[yftb$flag47b==1 & (! is.na(yftb$flag47b))]

i=0; 
for( myorf in as.character( orfs47b ) ) {
#for( myorf in as.character( orfs47b[1:5] ) ) {
 i = i + 1;
 cmd = paste ( 'touch /tmp/', i, '.txt', sep=''); #tracking the progress
 system( cmd );

 wt = c(); mut = c(); wt2 = c(); mut2 = c();

 tmp = w3.all[as.character(w3.all$orf) == myorf , ]
 tmp[1,c("orf","chr","gene")]
 out$chr[i]  = as.character( tmp$chr[1]);
 out$gene[i] = as.character( tmp$gene[1]);
 out$geneflag[i] = as.character( tmp$geneflag[1]);
 out$orf[i] = myorf;
 
 v1b = w1.all$FGMean2[as.character(w1.all$orf) == myorf ]; 
 v2b = w2.all$FGMean2[as.character(w2.all$orf) == myorf ]; 
 v3b = w3.all$FGMean2[as.character(w3.all$orf) == myorf ]; 
 wt2 = c(v1b,v2b,v3b);

 u1b = m1.all$FGMean2[as.character(m1.all$orf) == myorf ]; 
 u2b = m2.all$FGMean2[as.character(m2.all$orf) == myorf ]; 
 u3b = m3.all$FGMean2[as.character(m3.all$orf) == myorf ]; 
 mut2 = c(u1b,u2b,u3b);

 ##022008 change
 wt2 = ifelse( wt2<=0, 1, wt2);
 mut2 = ifelse( mut2<=0, 1, mut2);

 if ( max(c(wt2, mut2)> 1) )  { #avoid 1, 1, 1, situation
 t.comp = t.test( log2(wt2), log2(mut2) );  ## complete compensation
 t.1.5  = t.test( log2(wt2*1.5), log2(mut2), paired=T ); #zero compensation 
 yftb$p.complete.compensation[ yftb$ORF== orfs47b[i] ] = t.comp$p.value; 
 yftb$p.zero.compensation[ yftb$ORF== orfs47b[i] ] = t.1.5$p.value;  
 }
 
 my.ratio = yftb[  yftb$ORF== orfs47b[i], 7 ] 
 if ( t.comp$p.value >0.05 & t.1.5$p.value < 0.05 ) {
 #if ( t.comp$p.value >0.05 & t.1.5$p.value < 0.05 & my.ratio > 0.75 & my.ratio < 1.25 ) {
 	 yftb$compflag[  yftb$ORF== orfs47b[i] ] = 'compensated'; 
 	} else {
 	 yftb$compflag[  yftb$ORF== orfs47b[i] ] = 'not-compensated';
 	}
 
}


orfs5   = yftb$ORF[yftb$Ch==5 ]
i=0; 
for( myorf in as.character( orfs5 ) ) {
 i = i + 1;
 cmd = paste ( 'touch /tmp/', i, '.txt', sep=''); #tracking the progress
 system( cmd );

 wt = c(); mut = c(); wt2 = c(); mut2 = c();

 tmp = w3.all[as.character(w3.all$orf) == myorf , ]
 tmp[1,c("orf","chr","gene")]
 out$chr[i]  = as.character( tmp$chr[1]);
 out$gene[i] = as.character( tmp$gene[1]);
 out$geneflag[i] = as.character( tmp$geneflag[1]);
 out$orf[i] = myorf;
 
 v1b = w1.all$FGMean2[as.character(w1.all$orf) == myorf ]; 
 v2b = w2.all$FGMean2[as.character(w2.all$orf) == myorf ]; 
 v3b = w3.all$FGMean2[as.character(w3.all$orf) == myorf ]; 
 wt2 = c(v1b,v2b,v3b);

 u1b = m1.all$FGMean2[as.character(m1.all$orf) == myorf ]; 
 u2b = m2.all$FGMean2[as.character(m2.all$orf) == myorf ]; 
 u3b = m3.all$FGMean2[as.character(m3.all$orf) == myorf ]; 
 mut2 = c(u1b,u2b,u3b);

 ##022008 change
 wt2 = ifelse( wt2<=0, 1, wt2);
 mut2 = ifelse( mut2<=0, 1, mut2);

 if ( max(c(wt2, mut2)> 1) )  {
 t.comp = t.test( log2(wt2), log2(mut2) );  ## complete compensation
 t.1.5  = t.test( log2(wt2*0.5), log2(mut2), paired=T ); #zero compensation 
 yftb$p.complete.compensation[ yftb$ORF== orfs5[i] ] = t.comp$p.value; 
 yftb$p.zero.compensation[ yftb$ORF== orfs5[i] ] = t.1.5$p.value;  
 }

 my.ratio = yftb[  yftb$ORF== orfs5[i], 7 ] 
 #if ( t.comp$p.value >0.05 & t.1.5$p.value < 0.05 & my.ratio > 0.75 & my.ratio < 1.25 ) {
 if ( t.comp$p.value >0.05 & t.1.5$p.value < 0.05 ) {
 	 yftb$compflag[  yftb$ORF== orfs5[i] ] = 'compensated'; 
 	} else {
 	 yftb$compflag[  yftb$ORF== orfs5[i] ] = 'not-compensated';
 	}
}

table(yftb$compflag)

tmp = yftb[yftb$compflag=='compensated' & ( ! is.na(yftb$compflag)), ]
hist( tmp[,7])

tmp5 = yftb[yftb$Ch=='5', ]
tmp47b= yftb[yftb$flag47b==1 & ( ! is.na(yftb$compflag)), ]



write.csv(yftb, "_compensation.by.probe.050911b.csv")



#q("yes");


