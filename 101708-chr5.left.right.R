### Oct17,check 2fold down genes on left and right arms on chr5
## Oct10,08 after fixing chr annotaion erros
### based 0422308.chr.R
### using background correction and non-chr5 nomalization

rm(list=ls())

#  Exactly form of mulinomial sampling 
# vector of x and p
loglh.multinomial.sampling <- function( x, p ) {
  total = sum( x );
  if ( length(x) == length(p) ) { 
        y = p ^ x;
        ret <- lfactorial(total) - sum( lfactorial(x) ) + sum(log(y) );
  } else {
        ret <- NA;
  }
}

################################
#source ("source/021908-sor55-norm-merge.b.R");
# or 
# load ("021908.bkgrdCorrected.nonChr5Normalized.RData");
# load ("022008.bkgrdCorrected.nonChr5Normalized.RData");
 load( "101008.bkgrdCorrected.nonChr5Normalized.myChr.RData");

################# the expression results
cen.orfs = c('CEN1','CENR','CEN2','CEN3','CEN4','CEN5','CEN6','CEN7');
chrs     = c('Chr1','ChrR','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7');
locs	 = c( '1', 'R','2','3','4','5','6','7');

#tb2 = read.table("cal21-annotation.csv", sep="\t", header=T, as.is=TRUE, fill=T,
# colClasses = c("character","character","character","character","integer","integer","integer",
#"character") );

tb = read.delim("cal21-annotation.csv",
 colClasses = c("character","character","character","character","integer","integer","integer",
"character") );

x = tb$start < tb$end
head( cbind( x, as.character( tb$strand) ) ) # visual check, OK

cens = tb[ match(cen.orfs, tb$orf), ]; 

xscale = 3.2E6;  
xlim=c(0, xscale); 
yscale = 4;
ylim = c(0.1, 6);
yy = c(0.1, 0.2, 0.5, 1, 2, 4);
xx = pretty( range( tb$end, na.rm=T ) )

#cutoffs = quantile( out$fold.mut.by.wt, probs=c(0.05,0.95) );  ### 101308 change!!!!
cutoffs = c(0.6, 1.5); 

#window.bp = 0.5E5;  window = 17;  pcutoff = 1E-6;#too littel blue
#window.bp = 0.5E5;  window = 15;  pcutoff = 1E-6;
#window.bp = 0.5E5;  window = 25;  pcutoff = 1E-4;  #this too loose? 
#window.bp = 0.5E5;  window = 19;  pcutoff = 1E-4;  #this too loose? 

#window.bp = 0.5E5;  window = 15;  pcutoff = 1E-4;  #1 cluster for 2fold genes
window.bp = 0.5E5;  window = 15;  pcutoff = 1E-3;  # for 

expecteds = c(0.05, 0.9, 0.05) * window;

tmpheader = c("chr", "len", "clusLen");
fragmentout = data.frame( matrix(nrow= length(chrs), ncol=length(tmpheader) ) );
names(fragmentout) = tmpheader;

fragmentout$chr = chrs; 

  i=6
#for ( i in 1:length(chrs)) {
#for ( i in c(5,6,8) ) {

  	mychr = chrs[ i ];
	exp.chrtb = out[ grep( mychr, out$chr),] ;  ###change 100208, not work for >10 chrs
	
	cen.loc = cens$start[i]/2  + cens$end[i]/2
	exp.chrtb$LR = ifelse ( exp.chrtb$loc > cen.loc, 'R','L' );

	u = 1; 	#u =  mean( exp.chrtb$fold.mut.by.wt); # [1] 0.7891938
	
	for( j in 1:length(exp.chrtb$orf) ) {
		yin = exp.chrtb$fold.mut.by.wt[j];
		mycol = 'NA';
		if (( yin > cutoffs[2] ) && ( exp.chrtb$WtMean>5)) { mycol = 'blue';   ##
		} else if ( ( yin < cutoffs[1] ) & (yin > 0.4) ) { mycol = 'red';   ## change 101308
		} else { mycol= 'green'; }
		
		exp.chrtb$type[j] = mycol;
		y = ifelse ( yin > u, yin-u, yin-u );
	}
	
	sorted.locs = sort( exp.chrtb$loc );
	sliding.test = exp.chrtb[ match(sorted.locs, exp.chrtb$loc), ];
	sliding.test$num = 1:length(sliding.test[,1] );

	sliding.test[grep( "MTL", sliding.test$gene), ]
	sliding.test[ sliding.test$num<168 & sliding.test$num > 147, ]

	red.orfs = sliding.test$orf[ sliding.test$type == "red" ] #41 red orfs
	red.tb = sliding.test[ sliding.test$type == "red", ] #41 red orfs
	row.names(red.tb) = as.character(red.tb$num);
	pos = red.tb$num
	names(pos) = as.character( pos );
	
	x = table(exp.chrtb$LR);
	y = table(red.tb$LR);

	
	threshhold = 3; 
	hc = hclust( dist(pos), method="single" ) 
	clus.obs = cutree( hc, h= threshhold );
	obs = mean( table(clus.obs) );
	sort(table(clus.obs));
	
	#now permute, hclust, clus, 
	random = sliding.test;
	random$num = sample(random$num);
	redpos = random$num[random$type=="red"]
	hc = hclust( dist( redpos ), method="single");
	clus = cutree( hc, h=threshhold );
	random.clus.size = mean( table(clus) );

	single.run = function( threshhold ) {
		random$num = sample(random$num);
		redpos = random$num[random$type=="red"]
		hc = hclust( dist( redpos ), method="single");
		clus = cutree( hc, h=threshhold );
		random.clus.size = mean( table(clus) );
	}

	multi.run = function( n, threshhold ) {
	  x = c();
	  for ( i in 1:n ) {
	    x = c(x, single.run(threshhold)); 
	  }
	}	
	
	y = multi.run( 100, threshhold);
	
	hist(y);
	#Damn, this is not significant. 
			
quit("no");

	#sliding chi-sq test
	sliding.test = exp.chrtb[ match(sorted.locs, exp.chrtb$loc), ];
	sliding.test$clustertype = "green";            ### default is always green, Oct 3,08
	sliding.test$pchisq = NA;
	
	head(sliding.test);
	exp.chrtb[grep("orf19.362",exp.chrtb$orf),] #check passed
	
	for( j in 1:(length(sliding.test$loc) - window + 1) ) {
		xs = j:(j+window-1);
		mylocs = sliding.test$loc[xs];  mylocs = mylocs[ ! is.na(mylocs) ];
		
		if (( max(mylocs) - min(mylocs)) < window.bp ) {
			types = sliding.test$type[xs];
			downs = length( types[ types=="red"] );
			means = length( types[ types=="green"] );
			ups = length( types[ types=="blue"] );
			obs = c( downs, means, ups);
			obs[is.na(obs)] = 0;
			#expecteds = c(0.5, 9, 0.5);
			k2 = sum( (obs - expecteds)^2 / expecteds );
			center = floor( j + window/2 ); 
			sliding.test$pchisq[center] = pchisq( k2, length(obs)-1, lower.tail = F);
	
			if ( sliding.test$pchisq[ center ] < pcutoff ) {
				if ((obs[1]> expecteds[1]) && (obs[3]< expecteds[3]) ){
					sliding.test$clustertype[center] = "red"; 
				} else if ((obs[1] < expecteds[1]) && (obs[3] > expecteds[3]) ) {
					sliding.test$clustertype[center] = "blue"; 
				} 
			} else {
				sliding.test$clustertype[center] = "green"; #greens are not shown
			}
		
		} ## window.bp if loop
 
	}
	
	sliding.test[1:5, c(1,12:18)]
	tt = table( sliding.test$clustertype); tt;

	sliding.test[ sliding.test$clustertype == "red", ]
	
	
        #
	### calculate the block length
	myblockIndice = c(); 
	
	pt = 1; #current pointer
	while ( pt <= length(sliding.test[,1]) ) {
	      current = sliding.test$clustertype[pt];
	      #previous = sliding.test$clustertype[pt-1];
	      if ( ( current == "red" ) || (current=="blue") ) { 
	          myblockIndice = c( myblockIndice, pt );
	      }
 	      pt = pt + 1;
	}
	
	names(myblockIndice) = myblockIndice;
	hc = hclust( dist(myblockIndice), method="single" ) 
	#clus = cutree( hc, h=1);  ##over estimate the cluster numbers
	#clus = cutree( hc, h= floor(window/2));
	clus = cutree( hc, h= window);
	
	#sliding.test[95:111,c(1,3,7:9,15,17:20)] #1st cluster
	sliding.test[95:111,c(1,8:9,15,17:20)] #1st cluster, 5 reds

	sliding.test[213:239,c(1,8:9,15,17:20)] #2nd cluster 5 reds
	
	clustb = data.frame( cbind( as.integer(clus), as.integer(names(clus)) ) );
	clusLen = c();
	names( clustb ) = c("clus","index");
	for( nc in 1:length(unique(clus)) ) {
	   currentclus = clustb[clustb$clus==nc,]
	   start = min( currentclus$index );
	   end = max( currentclus$index );
	   adjust = floor( window/2 ); #overestimate clusters size, such as i, j when j-i < 7. 
	   #adjust = 1;
	   clusLen = c(clusLen, sliding.test$loc[end+adjust] - sliding.test$loc[start-adjust] ); #potential bug at boundries
	}
	
	fragmentout$len[i] = max( sliding.test$loc );
   	fragmentout$clusLen[i] = sum(clusLen); #to


#}#chr loop

fragmentout$clusFrag = fragmentout$clusLen / fragmentout$len;

write.csv(fragmentout, "_clusterLength.chr.101008.csv")

q("yes");

#########################################

	#start = min(sliding.test$loc);
	#end = max(sliding.test$loc);
	start = NA;
	end = NA;
	previousType ="green";
	
	   if ( is.na(sliding.test$clustertype[pt]) ) {
	      if ( (previousType == "red" ) || (previousType=="green") ) { 
	          ### end an old block
	          myblocks = rbind( myblocks, c(start, end) );
	          ### do not start a new block ????
	          start = NA; end = NA;
	      }
	   } else {
		if  ( sliding.test$clustertype[ pt ] == "blue" ) {
		  if (previousType == "blue")  {
			### extend the current block
			end = sliding.test$loc[pt];
		  } else { # previousType is "green". I don't consider blue-red change. 
			### the current block has ended, start a new block
			myblocks = rbind( myblocks, c(start, end) );
			start = sliding.test$loc[pt];
			end   = sliding.test$loc[pt];
		  }
		} else if (sliding.test$clustertype[pt] == "red") {
		  if (previousType == "red" ) {
			end = sliding.test$loc[pt];
		  } else { 
			myblocks = rbind( myblocks, c(start, end) );
			start = sliding.test$loc[pt];
			end   = sliding.test$loc[pt];
		  }
		} else { ### green 
		   start = NA; end = NA;
		} 
	   }
	   previousType = sliding.test$clustertype[ pt ]; 
	   previousType = ifelse( is.na(previousType), "missing", previousType ); 
	   pt = pt +1; 
	}
