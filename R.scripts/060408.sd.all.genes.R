### March 24, 08  find northern controls
###

 load ("022008.bkgrdCorrected.nonChr5Normalized.RData");

orfs  = as.character( unique( w1.all$orf ));

labels = c("orf","gene","FGMean","FGMean2","geneflag")
wt  = data.frame( matrix(nrow=1, ncol=length( labels ) ) );
names(wt)  = labels ;
mut = wt; 


for( myorf in as.character( orfs ) ) {
 
 v1 = w1.all[as.character(w1.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; v1;
 v2 = w2.all[as.character(w2.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; v2;
 v3 = w3.all[as.character(w3.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; v3;

 u1 = m1.all[as.character(m1.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; u1;
 u2 = m2.all[as.character(m2.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; u2;
 u3 = m3.all[as.character(m3.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; u3;
 
 wt  = rbind( wt,  v1, v2, v3);
 mut = rbind( mut, u1, u2, u3); 	

}

wt   = wt[ ! is.na(wt$orf), ]
mut = mut[ ! is.na(mut$orf), ]

tmp.tb = data.frame( wt$orf );
tmp.tb$ratio.mw = mut$FGMean2 / wt$FGMean2; 

for( i in 1:length(out[,1] ) ) {
#for( i in 1:10 ) {
  myorf = as.character( out$orf[i] );

  ratios = tmp.tb$ratio.mw[ tmp.tb$wt.orf == myorf ]
  out$sd.ratio[i]   = sd( ratios );
 
 }

#some genes are only on the last 2 chips, so no variance. 

write.table( out, "_out.sor55.nonCh5norm.with.sd.060408.csv", row.names=F);

q("yes");
