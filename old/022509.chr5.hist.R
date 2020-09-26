### plot 
### 102808 remove some orfs around MTLA locus
###oct11, remove background genes (bottom 5% genes)
###oct10, after fixing chr annotation erros
### 042408 partition chr5 genes

rm(list=ls())
 load( "101008.bkgrdCorrected.nonChr5Normalized.myChr.RData");

chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr6','Chr4','Chr7','Chr5');

wt.q  =  quantile( out$WtMean2,  probs = c(0.05) );
mut.q =  quantile( out$MutMean2, probs = c(0.05) );

#cutoffs.fold = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );
#cutoffs.mut  = quantile( out$MutMean2, probs=c(0.25,0.75) );

orfs.to.remove = c(
"orf19.3200",   # MTLA2
"orf19.3201",	# MTLA1	
"orf19.3197",	# PAP1	
"orf19.3199",	# PIKA	
"orf19.3198",	# OBPA	
"orf19.3982"	#orf19.3982
); #these are genes in MTL locus on Chr5 that Bulga wants to exclude

i = 8;
mychr = chrs[ i ];
tmp2 = out[out$chr== mychr,] ;
tmp3 = tmp2[ ! is.na(tmp2$orf), ]

tmp3$match = match( tmp3$orf, orfs.to.remove );
tmp4 = tmp3[ is.na(tmp3$match), ]#nice 501 genes

mtl = tmp3[ ! is.na(tmp3$match), ]
mtl = mtl[, -length(mtl[1,])]

# exp.chr5 = tmp[ ! is.na(tmp$orf), ] #507 genes
exp.chr5 = tmp4[, - length(tmp4[1,])]; #501 gens, 102808

tmp = exp.chr5[ exp.chr5$WtMean2 < wt.q[1] | exp.chr5$MutMean2 < mut.q[1] , ]
tmp = rbind( mtl, tmp);

write.csv( tmp, "_out.chr5.bkgrdGenes.022509.csv");

exp.chr5B = exp.chr5[ exp.chr5$WtMean2 > wt.q[1] & exp.chr5$MutMean2 > mut.q[1] , ]
#456 genes
write.csv(exp.chr5B, "_out.chr5.456genesAboveBkgrd.022509.csv");

left = seq( 0.1, 1.7, by=0.1); #left boundries
right = left + 0.1;  #right boundries
left[1] = -999999;
right[length(right)] = 9999;
out = data.frame( cbind(left, right) );

#exp.chr5 = exp.chr5B; 

for( i in 1:length(left) ) {
 sub = exp.chr5[ (exp.chr5$fold.mut.by.wt < right[i])&(exp.chr5$fold.mut.by.wt >= left[i] ),  ]; 
 out$freq[i] = length(sub[,1])
 
 subB = exp.chr5B[ (exp.chr5B$fold.mut.by.wt < right[i])&(exp.chr5B$fold.mut.by.wt >= left[i] ),  ]; 
 out$freqB[i] = length(subB[,1])
 
}

out$pos = out$left /2 + out$right /2 ;
out$pos[1] = 0.05;
out$pos[length(right)]=1.75;

pdf("_chr5.partition.020409.pdf", width=7,height=7);

plot( out$freq ~ out$pos, type='l',col="black", xlab='Ratio Sor55/3153A',ylab='Number of genes',xlim=c(0,4.5), las=2, axes=F );
points( out$freq ~ out$pos, pch=16);

cag1 = tmp3[tmp3$orf=="orf19.4015",]
x0 = cag1$fold.mut.by.wt[1];
y0 = 2 
arrows(x0,y0,x0,y0+50,code=1,lwd=1, col="red",length=0.1)
text(x0,y0+55,"CAG1");

axis(1,at=0:5);
axis(2,at=seq(0,100,by=20), las=2);

dev.off();

#lines( out$freq ~ out$pos, lty=2, col="blue" );
#points( out$freqB ~ out$pos, col="red" );
#lines( out$freqB ~ out$pos, lty=2, col="red" );
#legend( 1.1, 95, legend=c('all', 'w/o low expressed'), col=c('blue','red'), lty=c(2,2))


quit( "no" );
	
