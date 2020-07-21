### Dec 18 08  find semi-quantitative PCR controls
###

 load( "101008.bkgrdCorrected.nonChr5Normalized.myChr.RData");

orfs  = as.character( unique( w1.all$orf ));

labels = c("orf","gene","FGMean","FGMean2","geneflag")
#wt  = data.frame( matrix(nrow=1, ncol=length( labels ) ) );
#names(wt)  = labels ;
#mut = wt; 

#for( i in 1:10 ) {
for( i in 1:length(out$orf) ) {
 myorf = as.character( out$orf[i] );
 v1 = w1.all[as.character(w1.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; v1;
 v2 = w2.all[as.character(w2.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; v2;
 v3 = w3.all[as.character(w3.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; v3;

 u1 = m1.all[as.character(m1.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; u1;
 u2 = m2.all[as.character(m2.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; u2;
 u3 = m3.all[as.character(m3.all$orf) == myorf, c("orf","gene","FGMean","FGMean2","geneflag")]; u3;
 
 #wt  = rbind( wt,  v1, v2, v3);
 #mut = rbind( mut, u1, u2, u3); 	
 wt  <-  c(v1$FGMean2,v2$FGMean2,v3$FGMean2);
 mut <-  c(u1$FGMean2,u2$FGMean2,u3$FGMean2);
 wt[wt<1] = NA; 
 mut[mut<1]=NA;
 out$sd.wt[i]  = sd(wt, na.rm=T);
 out$sd.mut[i] = sd(mut, na.rm=T);
 out$sd.mut.vs.wt.ratio[i] = sd( mut/wt, na.rm=T );
}

write.csv( out, "_out.sor55.nonCh5norm.with.sd.121808.csv", row.names=F);
save( out, file="_121808.out.with.sd.RData");
load( "_121808.out.with.sd.RData" );

out = out[ ! is.na(out$sd.mut.vs.wt.ratio), ]

 
control = out[out$WtMean2<1500 & out$WtMean2>900 & out$sd.mut.vs.wt.ratio<0.2 & out$fold.mut.by.wt > 0.9 & out$fold.mut.by.wt<1.1, ] #36 genes

 write.csv(control,"_121808.control.for.ASU51.CSU53.csv");
 
control2 = out[out$WtMean2<380 & out$WtMean2>150 & out$sd.mut.vs.wt.ratio<0.5 & out$fold.mut.by.wt > 0.9 & out$fold.mut.by.wt<1.1, ] #only 2 genes
control2 = out[out$WtMean2<380 & out$WtMean2>150 & out$sd.mut.vs.wt.ratio<0.5, ]

 write.csv(control2,"_121808.control.for.ASU53.csv");


#some genes are only on the last 2 chips, so no variance. 

q("no");
