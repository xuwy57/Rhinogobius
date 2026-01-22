library(stringr)

bGrp <- T;
arrGrp1 <- c();
arrGrp2 <- c();

stdin<-file('stdin', 'r')

while(T) {
	sLn <- readLines(stdin, n=1);
	if (sLn == "END") {
		break;
	}

	arrNums <- as.numeric(unlist(str_split(sLn, "\t")));
	if (bGrp) {
		arrGrp1 <- arrNums;
	} else {
		arrGrp2 <- arrNums;
		cat(paste(mean(arrGrp1), mean(arrGrp2), median(arrGrp1), median(arrGrp2), wilcox.test(arrGrp1, arrGrp2, paired=F)$p.value, sep="\t"),"\n");
	}

	bGrp <- !bGrp;
}
