### Feb 11, 2008
### after background correction!!!
### histogram of expression ratio along chromosome, then histogram
### modified from 020508 version

rm(list=ls())

################################

#source ("source/021108-sor55-norm-merge.b.R");
# or 
 load ("021108.bkgrd.corrected.RData");


################# 
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7');

fields = c(  "chr", "WtMean2","MutMean2", "Mean.fold.mut.by.wt","WtMedian2","MutMedian2","Median.fold.mut.by.wt", "p");
out = data.frame( matrix(nrow=length(chrs), ncol= length(fields) ) );
names(out) = fields;

### average and medians by chromosome
for( c in 1:8) { #loop over 8 chromosomes
 chr = chrs[c]; 
 w1a.sub = w1.all$FGMean2[ as.character(w1.all$chr) == chr ];
 w2a.sub = w2.all$FGMean2[ as.character(w2.all$chr) == chr ];
 w3a.sub = w3.all$FGMean2[ as.character(w3.all$chr) == chr ];

 m1a.sub = m1.all$FGMean2[ as.character(m1.all$chr) == chr ];
 m2a.sub = m2.all$FGMean2[ as.character(m2.all$chr) == chr ];
 m3a.sub = m3.all$FGMean2[ as.character(m3.all$chr) == chr ];

 wt  = c( w1a.sub, w2a.sub, w3a.sub);
 mut = c( m1a.sub, m2a.sub, m3a.sub);

 out$chr[c] = chr;
 out$WtMean2[c] = mean( wt, na.rm=T );
 out$MutMean2[c] = mean( mut, na.rm=T );
 out$Mean.fold.mut.by.wt [c] = out$MutMean2[c] / out$WtMean2[c] 
 
 out$WtMedian2[c] = median( wt, na.rm=T );
 out$MutMedian2[c] = median( mut, na.rm=T );
 out$Median.fold.mut.by.wt [c] = out$MutMedian2[c] / out$WtMedian2[c] 

 t = t.test( wt, mut, paired=T);
 out$p[c] = t$p.value;

}


write.table(out, "_out.by.chr.021408.csv", sep="\t", row.name=F, quote=F);

q("no");



