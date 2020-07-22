### March 23, 08  find northern controls
###

 load ("022008.bkgrdCorrected.nonChr5Normalized.RData");

orfs  = c( "orf19.5007",
"orf19.5383",
"orf19.5928",
"orf19.5893",
"orf19.3651",
"orf19.3838",
"orf19.6814",
"orf19.4152",
"orf19.6745");

#labels = c( "orf", "chr", "gene", "geneflag", "WtMean","MutMean","WtMean2","MutMean2",  "fold.mut.by.wt", "p", "p.paired");
#out = data.frame( matrix(nrow=length( union(orfs, orfs2) ), ncol=length(labels) ) );
#names(out) = labels;

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


ctrl = data.frame( wt$orf );

ctrl$ratio.mw = mut$FGMean2 / wt$FGMean2; 

out2 = data.frame( orfs );
i = 1;
for( myorf in as.character(orfs) ) {

  genes = wt$gene[ wt$orf == myorf ]
  out2$gene[i] = genes[1];

  ratios = ctrl$ratio.mw[ ctrl$wt.orf == myorf ]
  out2$mean[i] = mean( ratios );
  out2$media[i] = median( ratios );
  out2$sd[i]   = sd( ratios );

  i = i + 1;
}

write.table( wt,   "_northern.control.wt.row.032308.csv", row.names=F);
write.table( mut, "_northern.control.mut.row.032308.csv", row.names=F);
write.table( out2, "_northern.controls.032308.csv", row.names=F);

q("no");
