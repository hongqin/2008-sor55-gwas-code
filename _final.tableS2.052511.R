### generate finale S2 table

rm(list=ls())

load( "030411.myannotation.bkgrdCorrected.nonChr457Normalized.RData" );
s2 = read.csv("s2.csv",header=T)
s2$ORF = as.character(s2$ORF)

pos= match( s2$ORF, out$orf)
s2$sd   = out$sd.wt[pos];
s2$sd.1 = out$sd.mut[pos];

s2[1:6,8]
s2[ s2$Mean<100 & s2$X3153A==1.0,8] = NA; 

write.csv(s2, "s2b.csv") 
