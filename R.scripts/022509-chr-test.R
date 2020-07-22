### Feb 21, 2008, add gene freq
### Feb 20, 2008 t-test by chr using log2 and histogram of ratio
### after background correction and non-chr5 nomalization
### histogram of expression ratio along chromosome, then histogram
### modified from 020508 version

rm(list=ls())

#source ("source/021908-sor55-norm-merge.b.R");
# or 
# load ("022008.bkgrdCorrected.nonChr5Normalized.RData");
 load ("101008.bkgrdCorrected.nonChr5Normalized.myChr.RData");

#remomve 5% as background genes for t-test and mean, median calculations
 wt.q  =  quantile( out$WtMean2,  probs = c(0.05) );
 mut.q =  quantile( out$MutMean2, probs = c(0.05) );

#################
cen.orfs = c('CEN1','CENR','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7');
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7');
locs	 = c( '1', 'R','2','3','4','5','6','7');

mat = matrix( seq(1,length(chrs)), nrow=length( chrs)/2, ncol= 2 ); 
layout(mat);

fields = c('chr', 'meanfold.sor55.by.wt', 'medianfold.sor55.by.wt', 'p.paired', 'gene.num');
out.chr = data.frame( matrix(nrow=length(chrs), ncol= length(fields) ) );
names( out.chr ) = fields; 

#i=4;
for ( i in 1:length(chrs)) {
	mychr = chrs[ i ];
	tmp = out[out$chr== mychr,] ;
	tmp = tmp[ ! is.na(tmp[,1]), ]

	out.chr$gene.num[i] = length( tmp$chr );
	
	#for t-test 
        tmp$WtMean2   = ifelse( tmp$WtMean2<= wt.q[1],  NA, tmp$WtMean2);  #remove bottom 5% genes
        tmp$MuttMean2 = ifelse( tmp$MutMean2<=mut.q[1], NA, tmp$MutMean2); #remove bottom 5% genes

	t = t.test(log2(tmp$WtMean2), log2(tmp$MutMean2), paired=T);
	out.chr$p.paired[i] =  t$p.value; 
	
	# for histogram
	tmp$fold.mut.by.wt[tmp$fold.mut.by.wt==1]= NA;
	exp.chrtb = tmp[ ! is.na(tmp$orf), ]
	exp.chrtb = exp.chrtb[ ! is.na(exp.chrtb$fold.mut.by.wt), ]
	
	out.chr$chr[i] = mychr;
	out.chr$meanfold.sor55.by.wt[i] = mean(exp.chrtb$fold.mut.by.wt);
	out.chr$medianfold.sor55.by.wt[i] = median(exp.chrtb$fold.mut.by.wt);
	
	par(mar=c(2,2,1,1));
	hist( log2(exp.chrtb$fold.mut.by.wt),col='blue',border='NA',br=50,xlim=c(-5,5),
	 ylab='',
	 main= paste( "log2(sor55/wt) at ", mychr, sep=' '));
	 	
		
}#chr loop
#dev.off();

out.chr;

write.table(out.chr, "_out.sor55.by.chr.022509.csv", sep="\t", row.name=F, quote=F);

#q("no");

