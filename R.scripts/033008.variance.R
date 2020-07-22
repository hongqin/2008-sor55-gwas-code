### March 30, 08  Bulgac insist upon all values, so I will email her all six arrays

 load ("022008.bkgrdCorrected.nonChr5Normalized.RData");

write.table( w1.all, "normalized.arrrays/_w1.csv", row.names=F, quote=F, sep="\t");
write.table( w2.all, "normalized.arrrays/_w2.csv", row.names=F, quote=F, sep="\t");
write.table( w3.all, "normalized.arrrays/_w3.csv", row.names=F, quote=F, sep="\t");

write.table( m1.all, "normalized.arrrays/_m1.csv", row.names=F, quote=F, sep="\t");
write.table( m2.all, "normalized.arrrays/_m2.csv", row.names=F, quote=F, sep="\t");
write.table( m3.all, "normalized.arrrays/_m3.csv", row.names=F, quote=F, sep="\t");


q("no");



