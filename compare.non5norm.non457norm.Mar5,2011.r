old = read.delim("_out.sor55.nonCh5norm.with.sd.060408.csv");
new = read.delim("_out.sor55.nonCh4,5,7norm.030411.csv");
old$orf = as.character(old$orf)
new$orf = as.character(new$orf)
rownames(old) = old$orf
rownames(new) = new$orf

myorfs = intersect( old$orf, new$orf )
plot( new[myorfs, 'fold.mut.by.wt'] ~ old[myorfs, 'fold.mut.by.wt'], xlim=c(0,4), ylim=c(0,4))
#generally agree with each other. 

#summary(lm(  ))

#plot(new$fold.mut.by.wt[myorfs] ~ old$fold.mut.by.wt[myorfs])

#_out.sor55.by.chr.031809.csv
#_out.sor55.nonCh5norm.083109.csv
#_out.sor55.nonCh5norm.092909.csv
#_out.sor55.nonCh5norm.with.sd.030311.csv
#_out.sor55.nonCh5norm.with.sd.121808.csv

all = merge(old, new)
