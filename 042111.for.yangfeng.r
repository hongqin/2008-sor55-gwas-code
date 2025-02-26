### script to generate normalized array for yang feng. April 21, 2011

### non4,5,7 normalization
### Sep 29, 09, change negative values to positives; 
### Aug 31-Sep1, 09 update annotation based on cal21-annotation.csv
### Feb 20, 08, log2 transformation for t.test
### Feb 19, 2008 norm by non-chr5 genes
### background correction by chips, merge results, calculate sor55/wt
### use non-chr5 orfs as normalization, redo with data with re-annotated CSU-ASU genes

rm(list=ls())
files = c(
"_reformated.3153A_4006873_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv",  #1 wt uneven
"_reformated.Results_SolX_CMI.txt.csv",	#2 wt
"_reformated.Sol125a_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #3 mut
"_reformated.Sol125b_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #4 mut
"_reformated.SolXb_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #5  wt
"_reformated.Sor55_4006777_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv", #6  mut uneven
"_reformated.3153A-121407.csv.csv", # 7 wt
"_reformated.Sor55-121407.csv.csv"  #8  mut
);

wildtype.files  = files[c(2, 5, 7, 1)];  # control files
mutant.files    = files[c(3, 4, 8, 6)];  # experiment files

### read in control data
w1.all = read.delim( wildtype.files[1], header=T, sep="\t", fill=T); 
w2.all = read.delim( wildtype.files[2], header=T, sep="\t", fill=T); 
w3.all = read.delim( wildtype.files[3], header=T, sep="\t", fill=T); 
#w4.all = read.delim( wildtype.files[4], header=T, sep="\t", fill=T); #need for CS

w1.all = w1.all[ ! is.na(w1.all$geneflag), ]
w2.all = w2.all[ ! is.na(w2.all$geneflag), ]
w3.all = w3.all[ ! is.na(w3.all$geneflag), ]
#w4.all = w4.all[ ! is.na(w4.all$geneflag), ]

#pw1.all = w1.all[ w1.all$geneflag==0, ]
#pw2.all = w2.all[ w2.all$geneflag==0, ]
#pw3.all = w3.all[ w3.all$geneflag==0, ]
#pw1.all = pw1.all[ ! is.na(pw1.all$geneflag), ]
#pw2.all = pw2.all[ ! is.na(pw2.all$geneflag), ]
#pw3.all = pw3.all[ ! is.na(pw3.all$geneflag), ]

### read in experimental data
m1.all = read.delim( mutant.files[1], header=T, sep="\t", fill=T); 
m2.all = read.delim( mutant.files[2], header=T, sep="\t", fill=T); 
m3.all = read.delim( mutant.files[3], header=T, sep="\t", fill=T); 
#m4.all = read.delim( mutant.files[4], header=T, sep="\t", fill=T); 

m1.all = m1.all[ ! is.na(m1.all$geneflag), ];
m2.all = m2.all[ ! is.na(m2.all$geneflag), ];
m3.all = m3.all[ ! is.na(m3.all$geneflag), ];
#m4.all = m3.all[ ! is.na(m4.all$geneflag), ];

#pm1.all = m1.all[ m1.all$geneflag==0, ];
#pm2.all = m2.all[ m2.all$geneflag==0, ];
#pm3.all = m3.all[ m3.all$geneflag==0, ];
#pm1.all = pm1.all[ ! is.na(pm1.all$geneflag), ];
#pm2.all = pm2.all[ ! is.na(pm2.all$geneflag), ];
#pm3.all = pm3.all[ ! is.na(pm3.all$geneflag), ];

### update the annotations ### 083109 change
cen.orfs = c('CEN1','CENR','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7');
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7');
locs     = c( '1', 'R','2','3','4','5','6','7');

dic= data.frame( cbind( chrs, locs) )
dic$chrs = as.character( dic$chrs );
dic$locs = as.character( dic$locs );

#
tb = read.delim("cal21-annotation.csv");
antb = read.delim("cal21-annotation.csv", colClasses = c("character","character","character","character","integer","integer","integer", "character") );

#tb$chr = as.character( tb$chr );
cens = antb[ match(cen.orfs, antb$orf), ]; 

# combine table 
update.table = function( intb ) {
 intb$chr.combimatrix = intb$chr; #store the original chr annotation
 intb$start = antb$start[ match( intb$orf, antb$orf ) ]; # 100208
 intb$end   = antb$end[ match( intb$orf, antb$orf ) ]; #  100208       
 intb$strand= antb$strand[ match( intb$orf, antb$orf ) ]; #  100208    
 intb$chr2  = antb$chr[ match( intb$orf, antb$orf ) ]; #  100208       
 intb$chr   = dic$chrs[ match( intb$chr2, dic$locs) ] 
 intb$loc   = intb$start /2  + intb$end/2 
 ret <- intb;
}

w1.all = update.table(w1.all);
w2.all = update.table(w2.all);
w3.all = update.table(w3.all);
#w4.all = update.table(w4.all);
m1.all = update.table(m1.all);
m2.all = update.table(m2.all);
m3.all = update.table(m3.all);
#m4.all = update.table(m4.all);

several = ceiling( runif(5)*length(w1.all[,1]) )
head(w1.all[ several, c("chr","chr2", "chr.combimatrix")])
head(m1.all[ several, c("chr","chr2", "chr.combimatrix")])
head(m3.all[ several, c("chr","chr2", "chr.combimatrix")])

### background corection
w1.background   = median( w1.all$FGMean[w1.all$geneflag==0])
w2.background   = median( w2.all$FGMean[w2.all$geneflag==0])
w3.background   = median( w3.all$FGMean[w3.all$geneflag==0])
#w4.background   = median( w4.all$FGMean[w4.all$geneflag==0])
#c(w1.background,w2.background,w3.background ,w4.background )


m1.background   = median( m1.all$FGMean[m1.all$geneflag==0])
m2.background   = median( m2.all$FGMean[m2.all$geneflag==0])
m3.background   = median( m3.all$FGMean[m3.all$geneflag==0])
#m4.background   = median( m4.all$FGMean[m4.all$geneflag==0])
#c(m1.background,m2.background,m3.background ,m4.background )

w1.all$FGMean2 = w1.all$FGMean - w1.background;
w2.all$FGMean2 = w2.all$FGMean - w2.background;
w3.all$FGMean2 = w3.all$FGMean - w3.background;
#w4.all$FGMean2 = w4.all$FGMean - w4.background;
w1.all$FGMean2[w1.all$FGMean2<0] = 0;  #092909 change
w2.all$FGMean2[w2.all$FGMean2<0] = 0; 
w3.all$FGMean2[w3.all$FGMean2<0] = 0; 

m1.all$FGMean2 = m1.all$FGMean - m1.background;
m2.all$FGMean2 = m2.all$FGMean - m2.background;
m3.all$FGMean2 = m3.all$FGMean - m3.background;
#m4.all$FGMean2 = m4.all$FGMean - m4.background;
m1.all$FGMean2[m1.all$FGMean2<0] = 0; 
m2.all$FGMean2[m2.all$FGMean2<0] = 0; 
m3.all$FGMean2[m3.all$FGMean2<0] = 0; 

### normalization, non4,5,7
#scale = sum( w1.all$FGMean2[w1.all$geneflag==1], na.rm=T);
scale.5 = sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)], na.rm=T); #030411 change
scale = sum( w1.all$FGMean2[(w1.all$chr != 'Chr4') & (w1.all$chr != 'Chr5') & (w1.all$chr != 'Chr7') & (w1.all$geneflag==1)], na.rm=T); #030411 change

w2.all$FGMean2 = w2.all$FGMean2 * scale / sum( w2.all$FGMean2[(w2.all$chr != 'Chr4') & (w2.all$chr != 'Chr5') & (w2.all$chr != 'Chr7') & (w2.all$geneflag==1)],na.rm=T);

w3.all$FGMean2 = w3.all$FGMean2 * scale / sum( w3.all$FGMean2[(w3.all$chr != 'Chr4') & (w3.all$chr != 'Chr5') & (w3.all$chr != 'Chr7') & (w3.all$geneflag==1)],na.rm=T);

#w4.all$FGMean3 = w4.all$FGMean * scale / sum( w4.all$FGMean2[(w4.all$chr != 'Chr5') & (w4.all$geneflag==1)],na.rm=T);

m1.all$FGMean2 = m1.all$FGMean2 * scale / sum( m1.all$FGMean2[(m1.all$chr != 'Chr4') & (m1.all$chr != 'Chr5') & (m1.all$chr != 'Chr7') & (m1.all$geneflag==1)],na.rm=T);

m2.all$FGMean2 = m2.all$FGMean2 * scale / sum( m2.all$FGMean2[(m2.all$chr != 'Chr4') & (m2.all$chr != 'Chr5') & (m2.all$chr != 'Chr7') & (m2.all$geneflag==1)],na.rm=T);

m3.all$FGMean2 = m3.all$FGMean2 * scale / sum( m3.all$FGMean2[(m3.all$chr != 'Chr4') & (m3.all$chr != 'Chr5') & (m3.all$chr != 'Chr7') & (m3.all$geneflag==1)],na.rm=T);

#m4.all$FGMean2 = m4.all$FGMean * scale / sum( m4.all$FGMean2[(m4.all$chr != 'Chr5') & (m4.all$geneflag==1)],na.rm=T);

### do some checks
w3.all[w3.all$orf=="orf19.2896",c("FGMean","FGMean2")] #SOU1 has 12 probes

w2.pos = (w2.all$chr != 'Chr4') & (w2.all$chr != 'Chr5') & (w2.all$chr != 'Chr7') & (w2.all$geneflag==1)
w3.pos = (w3.all$chr != 'Chr4') & (w3.all$chr != 'Chr5') & (w3.all$chr != 'Chr7') & (w3.all$geneflag==1)
m1.pos = (m1.all$chr != 'Chr4') & (m1.all$chr != 'Chr5') & (m1.all$chr != 'Chr7') & (m1.all$geneflag==1)
m2.pos = (m2.all$chr != 'Chr4') & (m2.all$chr != 'Chr5') & (m2.all$chr != 'Chr7') & (m2.all$geneflag==1)
m3.pos = (m3.all$chr != 'Chr4') & (m3.all$chr != 'Chr5') & (m3.all$chr != 'Chr7') & (m3.all$geneflag==1)

#check1 = sum( w2.all$FGMean2[(w2.all$chr != 'Chr5') & (w2.all$geneflag==1)],na.rm=T) / sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)],na.rm=T)
#check2 = sum( w3.all$FGMean2[(w3.all$chr != 'Chr5') & (w3.all$geneflag==1)],na.rm=T) / sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)],na.rm=T)
#checkm1 = sum( m1.all$FGMean2[(m1.all$chr != 'Chr5') & (m1.all$geneflag==1)],na.rm=T) / sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)],na.rm=T)
#checkm2 = sum( m2.all$FGMean2[(m2.all$chr != 'Chr5') & (m2.all$geneflag==1)],na.rm=T) / sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)],na.rm=T)
#checkm3 = sum( m3.all$FGMean2[(m3.all$chr != 'Chr5') & (m3.all$geneflag==1)],na.rm=T) / sum( w1.all$FGMean2[(w1.all$chr != 'Chr5') & (w1.all$geneflag==1)],na.rm=T)

check2 =sum(w2.all$FGMean2[w2.pos], na.rm=T) /scale 
check3 =sum(w3.all$FGMean2[w3.pos], na.rm=T) /scale 
checkm1 =sum(m1.all$FGMean2[m1.pos], na.rm=T) /scale 
checkm2 =sum(m2.all$FGMean2[m2.pos], na.rm=T) /scale 
checkm3 =sum(m3.all$FGMean2[m3.pos], na.rm=T) /scale 
c(check2, check3, checkm1, checkm2, checkm3); #these should be 1, 1, 1, 1, 1

write.csv(w1.all, 'for.yangfeng.042111/w1.csv')
write.csv(w2.all, 'for.yangfeng.042111/w2.csv')
write.csv(w3.all, 'for.yangfeng.042111/w3.csv')
write.csv(m1.all, 'for.yangfeng.042111/m1.csv')
write.csv(m2.all, 'for.yangfeng.042111/m2.csv')
write.csv(m3.all, 'for.yangfeng.042111/m3.csv')

quit('no')

########################################
### found out all genes 
orfs  = as.character( unique( w1.all$orf[ (! is.na(w1.all$geneflag)) & ( w1.all$geneflag == 1) ] ));
orfs = orfs[ -grep("wrong", orfs) ] #090109 remove the wrong orf

orfs2  = as.character( unique( w3.all$orf[ (! is.na(w3.all$geneflag)) & ( w3.all$geneflag == 1) ] ));
orfs.new = setdiff( orfs2, orfs);
#tmp = union( orfs, orfs2);
#write.table(orfs, "orfs.092909.csv", col.names=T, row.names=T, quote=F, sep="\t");
##  orfs = orfs[1:100]  #debug

labels = c( "orf", "chr", "gene", "geneflag", "WtMean","MutMean","WtMean2","MutMean2", "sd.wt","sd.mut","fold.mut.by.wt", "sd.ratio","p", "p.paired");
out = data.frame( matrix(nrow=length( orfs ), ncol=length(labels) ) );
#out = data.frame( matrix(nrow=length( union(orfs, orfs2) ), ncol=length(labels) ) );
names(out) = labels;

labels2 = c( "orf", "chr", "gene", "w1", "w2", "w3", "m1", "m2", "m3" );
data = data.frame( matrix( nrow=1, ncol=length(labels2) ) );
names(data) = labels2;

 i = 0;
for( myorf in as.character( orfs ) ) {
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
 
 v1 = w1.all$FGMean[as.character(w1.all$orf) == myorf ]; v1;
 v2 = w2.all$FGMean[as.character(w2.all$orf) == myorf ]; v2;
 v3 = w3.all$FGMean[as.character(w3.all$orf) == myorf ]; v3;
 wt = c(v1,v2,v3);

 u1 = m1.all$FGMean[as.character(m1.all$orf) == myorf ]; u1;
 u2 = m2.all$FGMean[as.character(m2.all$orf) == myorf ]; u2;
 u3 = m3.all$FGMean[as.character(m3.all$orf) == myorf ]; u3;
 mut = c(u1,u2,u3);

 out$WtMean[i]  = mean( wt );
 out$MutMean[i] = mean( mut );

 v1b = w1.all$FGMean2[as.character(w1.all$orf) == myorf ]; 
 v2b = w2.all$FGMean2[as.character(w2.all$orf) == myorf ]; 
 v3b = w3.all$FGMean2[as.character(w3.all$orf) == myorf ]; 
 wt2 = c(v1b,v2b,v3b);

 u1b = m1.all$FGMean2[as.character(m1.all$orf) == myorf ]; 
 u2b = m2.all$FGMean2[as.character(m2.all$orf) == myorf ]; 
 u3b = m3.all$FGMean2[as.character(m3.all$orf) == myorf ]; 
 mut2 = c(u1b,u2b,u3b);

 tmpw1 = w1.all[as.character(w1.all$orf) == myorf, c("orf", "chr", "gene", "FGMean2") ]
 tmpw2 = w2.all[as.character(w2.all$orf) == myorf, c( "FGMean2") ]
 tmpw3 = w3.all[as.character(w3.all$orf) == myorf, c( "FGMean2") ]
 tmpm1 = m1.all[as.character(m1.all$orf) == myorf, c( "FGMean2") ]
 tmpm2 = m2.all[as.character(m2.all$orf) == myorf, c( "FGMean2") ]
 tmpm3 = m3.all[as.character(m3.all$orf) == myorf, c( "FGMean2") ]
 
 rownames(tmpw1) = c();
 tmpdata = cbind( tmpw1, tmpw2, tmpw3, tmpm1, tmpm2, tmpm3 );
 names(tmpdata) = labels2; 
 data = rbind( data, tmpdata);

 #out[i,c("w1", "w2", "w3")] = wt2;  #090109 change, wrong
 #out[i,c("m1", "m2", "m3")] = mut2; #090109 change,wrong, due to multuple probes

 out$WtMean2[i]  = mean( wt2 );
 out$MutMean2[i] = mean( mut2 );

 ##022008 change
 wt2 = ifelse( wt2<=0, 1, wt2);
 mut2 = ifelse( mut2<=0, 1, mut2);

 t = t.test( log2(wt2), log2(mut2) );  ##022008 change
 t2 = t.test( log2(wt2), log2(mut2), paired=T);    ##022008 change
 out$p[i] = t$p.value; 
 out$p.paired[i] = t2$p.value;
 
 out$sd.wt[i]  = sd(wt2, na.rm=T);
 out$sd.mut[i] = sd(mut2, na.rm=T);
 out$sd.ratio[i] = sd( mut2/wt2, na.rm=T );
 
 #sd is calculated before truncation
 #sd.fold.square =  var(wt2)/ mean(wt2)^2 + var(mut2) / mean(mut2)^2 ;
 
 if ( mean(wt2) > 100 ) { #083109 this was later changed to 200 after averaging
     tmp = ifelse (mean(mut2)>1, mean(mut2), 1);
     out$fold.mut.by.wt[i] = tmp / mean( wt2 );    
 } else { #wt expression < 100
   if(  mean(mut2) > 100 ) {
     tmpwt = ifelse (mean(wt2)>1, mean(wt2), 1);
     out$fold.mut.by.wt[i] = mean(mut2) / tmpwt;   
   } else if ((mean(mut2)<10)&&(mean(wt2)>50) ) {
     tmp = ifelse (mean(mut2)>0, mean(mut2), 1);
     out$fold.mut.by.wt[i] = tmp / mean( wt2 );   
   } else { #wt and mut are weak
     out$fold.mut.by.wt[i] = 1 / 1;       #veiw as unchanged.
   } 
 }

 write.table(data, "_probe.data.merged.sor55.nonCh4,5,7norm.030411.csv", col.names=T, row.names=F, quote=F, sep="\t");
 #ch5.tb = data[ (data$chr=='Chr5') & (! is.na(data$chr)), ]
 #write.table(ch5.tb, "_chr5.probe.data.nonCh5norm.092909.csv", col.names=T, row.names=F, quote=F, sep="\t");
}

summary(out);

write.table(out, "_out.sor55.nonCh4,5,7norm.030411.csv", col.names=T, row.names=F, quote=F, sep="\t");

### 
 #up = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt > 1.3), ]
 #up = up[ up$geneflag==1,]
 #write.table(up, "_out.sor55.up.083109.csv", col.names=T, row.names=F, quote=F, sep="\t");

 #down = out[ (out$p.paired<0.01) & (out$fold.mut.by.wt < 0.5), ]
 #write.table(down, "_out.sor55.down.083109.csv", col.names=T, row.names=F, quote=F, sep="\t");

### bottom quatiles for "non-expressed" genes
 #out.qc = out[out$geneflag==0, ]

 #cutoff.wt  =  quantile( out$WtMean2,  probs=c(0.05,0.95), na.rm=T );
 #cutoff.mut =  quantile( out$MutMean2, probs=c(0.05,0.95), na.rm=T );

 #bottom = out[ ((out$WtMean2 < cutoff.wt[1]) & (out$MutMean2 < cutoff.mut[1])), ]
 #tmp = bottom[bottom$geneflag==0,]
 #tmp; #there are 25 control here, 

#write.table(bottom, "_out.sor55.bottom5percent.083109.csv", col.names=T, row.names=F, quote=F, sep="\t");

save.image( "030411.myannotation.bkgrdCorrected.nonChr457Normalized.RData" );

q("yes");


for( myorf in as.character( orfs.new ) ) {
 i = i + 1;
 wt = c(); mut = c(); wt2 = c(); mut2 = c();

 tmp = w3.all[as.character(w3.all$orf) == myorf , ]
 tmp[1,c("orf","chr","gene")]
 out$chr[i]  = as.character( tmp$chr[1]);
 out$gene[i] = as.character( tmp$gene[1]);
 out$geneflag[i] = as.character( tmp$geneflag[1]);
 out$orf[i] = myorf;
 
 v3 = w3.all$FGMean[as.character(w3.all$orf) == myorf ]; v3;
 v4 = w4.all$FGMean[as.character(w4.all$orf) == myorf ]; v3;
 wt = c(v3,v4);

 u3 = m3.all$FGMean[as.character(m3.all$orf) == myorf ]; u3;
 u4 = m4.all$FGMean[as.character(m4.all$orf) == myorf ]; u3;
 mut = c(u3,u4);

 out$WtMean[i]  = mean( wt );
 out$MutMean[i] = mean( mut );

 v3 = w3.all$FGMean2[as.character(w3.all$orf) == myorf ]; v3;
 v4 = w4.all$FGMean2[as.character(w4.all$orf) == myorf ]; v3;
 wt2 = c(v3,v4);

 u3 = m3.all$FGMean2[as.character(m3.all$orf) == myorf ]; u3;
 u4 = m4.all$FGMean2[as.character(m4.all$orf) == myorf ]; u3;
 mut2 = c(u3,u4);

 out$WtMean2[i]  = mean( wt2 );
 out$MutMean2[i] = mean( mut2 );
 
 ##022008 change
 wt2 = ifelse( wt2<=0, 1, wt2);
 mut2 = ifelse( mut2<=0, 1, mut2);

 t = t.test( log2(wt2), log2(mut2) );  ##022008 change
 t2 = t.test( log2(wt2), log2(mut2), paired=T);  ##022008 change
 out$p[i] = t$p.value; 
 out$p.paired[i] = t2$p.value;
 
 if ( mean(wt2) > 100 ) {
     tmp = ifelse (mean(mut2)>1, mean(mut2), 1);
     out$fold.mut.by.wt[i] = tmp / mean( wt2 );    
 } else { #wt expression < 100
   if(  mean(mut2) > 100 ) {
     tmpwt = ifelse (mean(wt2)>1, mean(wt2), 1);
     out$fold.mut.by.wt[i] = mean(mut2) / tmpwt;   
   } else if ((mean(mut2)<10)&&(mean(wt2)>50) ) {
     tmp = ifelse (mean(mut2)>0, mean(mut2), 1);
     out$fold.mut.by.wt[i] = tmp / mean( wt2 );   
   } else { #wt and mut are weak
     out$fold.mut.by.wt[i] = 1;       #veiw as unchanged.
   } 
 }

}

