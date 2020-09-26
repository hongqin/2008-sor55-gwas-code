s2 =read.csv( "s2.csv")
ch5 = s2[ s2$A21_CHROMOSOME=='Ca21chr5', ]

ch5$weak = 0;
ch5$weak[ch5$WtMean2 < 100 & ch5$MutMean2 < 100] = 1;

top = ch5[ ch5$weak==0, ]
bottom = ch5[ ch5$weak==1, ]

new = rbind( top, bottom)
write.csv(new, "reordered.ch5.csv,040611.csv")


