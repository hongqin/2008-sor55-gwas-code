### March 23, 2011 update chr location to output

### March 5, 2011 use non4,5,7 normalization

### 101008 use my annotation, there are errors in Combimatrix's chr annotation files, such as orf19.5320

### 042308 remove hist part. remove x-axis. tick inside
### based 021908-chr-hist.R
### after background correction and non-chr5 nomalization
### histogram of expression ratio along chromosome, then histogram
### modified from 020508 version

rm(list=ls())


################################
 load ("030411.myannotation.bkgrdCorrected.nonChr457Normalized.RData");

cen.orfs = c('CEN1','CENR','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7');
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7');
locs	 = c( '1', 'R','2','3','4','5','6','7');

dic= data.frame( cbind( chrs, locs) )
dic$chrs = as.character( dic$chrs );
dic$locs = as.character( dic$locs );

tb = read.delim("cal21-annotation.csv",
 colClasses = c("character","character","character","character","integer","integer","integer",
"character") );
cens = tb[ match(cen.orfs, tb$orf), ]; 

## now combine tables   101008 change
out$start = tb$start[ match( out$orf, tb$orf ) ]; # 100208
out$end   = tb$end[ match( out$orf, tb$orf ) ]; #  100208	
out$strand= tb$strand[ match( out$orf, tb$orf ) ]; #  100208	
out$chr2  = tb$chr[ match( out$orf, tb$orf ) ]; #  100208	
out$chr   = dic$chrs[ match( out$chr2, dic$locs) ] 
out$loc   = out$start /2  + out$end/2

table(out$chr); table( out$chr2); #check, OK

tb4 = out[out$chr=='Chr4',]
tb4 = tb4[ ! is.na(tb4$orf), ]
tb7 = out[out$chr=='Chr7',]
tb7 = tb7[ ! is.na(tb7$orf), ]

tb4 = tb4[order(tb4$loc), ]
tb7 = tb7[order(tb7$loc), ]

write.csv(out, "_out.sor55.nonCh4,5,7norm.032311.csv")
write.csv(tb4, "_out.sor55.nonCh4,5,7norm.Chr4.032311.csv")
write.csv(tb7, "_out.sor55.nonCh4,5,7norm.Chr7.032311.csv")


q("yes");

