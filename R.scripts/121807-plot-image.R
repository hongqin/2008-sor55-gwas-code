# Bulgac redo 2 hybridizations using new chips

rm(list=ls())

files = c(
"_reformated.3153A_4006873_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv",  #1 wt uneven
"_reformated.Results_SolX_CMI.txt.csv",		#2 wt
"_reformated.Sol125a_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #3 mut
"_reformated.Sol125b_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #4 mut
"_reformated.SolXb_Cyanine5_80focus_5um_40pmt_Feature_Data.txt.csv", #5  wt
"_reformated.Sor55_4006777_Cy5_80focus_90LP_50PMT_5micron_Feature_Data.txt.csv", #6  mut uneven
"_reformated.3153A-121407.csv.csv", # 7 wt
"_reformated.Sor55-121407.csv.csv"  #8  mut
);

#selected.files    = files[c(1,6, 2,3, 5,4, 7,8 )];  # all data
selected.files    = files[c(      2,3, 5,4, 7,8 )];  # only the good files

pdf ( paste("_image.sorbose.121807.pdf"), width=10 ) 

for( i in 1:length(selected.files) ) {

w1.all = read.delim( selected.files[i], header=T, sep="\t", fill=T); 

max(w1.all$Row, na.rm=T); #224 rows
max(w1.all$Column, na.rm=T); #56 columns
#unique(w1.all$orf);  #6372 levels

w1.mat = matrix( nrow=224, ncol=56);
for( r in 1:length(w1.all[,1] ) ) {
  if ( ( ! is.na(w1.all$Row[r]) ) & ( ! is.na(w1.all$Column[r]) ) ) {
   w1.mat[ w1.all$Row[r], w1.all$Column[r] ] = w1.all$FGMean[r];
  }
}

# w1.mat = ifelse( w1.mat>ceilings[i], ceilings[i], w1.mat);
image( 1:224, 1:56, w1.mat );
title ( selected.files[i]  )
}

dev.off();

q("no");
