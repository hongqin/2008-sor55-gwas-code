s2 =read.csv( "s2.csv")
ch6 = s2[ s2$A21_CHROMOSOME=='Ca21chr6', ]

ch6$weak = 0;
ch6$weak[ch6$WtMean2 < 100 & ch6$MutMean2 < 100] = 1;

top = ch6[ ch6$weak==0, ]
bottom = ch6[ ch6$weak==1, ]

new = rbind( top, bottom)

write.csv(new, "reordered.ch6.csv,041111.csv")


