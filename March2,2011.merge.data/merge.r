ch4 = read.csv("ch4.tb1.csv", header=F)
ch7 = read.csv("ch7.tb1.csv", header=F)
ch4[,1] = as.character(ch4[,1])
ch7[,1] = as.character(ch7[,1])

all = read.csv("all2.060408.csv")
all$orf = as.character(all$orf)

tmp4 = intersect( all$orf, ch4[,1]) #test
tmp7 = intersect( all$orf, ch7[,1]) #test

match4 = match(ch4[,1], all$orf) #match positions
match7 = match(ch7[,1], all$orf)

#merge
ch4b = data.frame( cbind(ch4, all2[match4,]) ) 
ch7b = data.frame( cbind(ch7, all2[match7,]) )

write.csv(ch4b, "ch4.March2,2011.csv");
write.csv(ch7b, "ch7.March2,2011.csv");

