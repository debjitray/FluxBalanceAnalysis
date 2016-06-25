library(Biobase)
library (Hmisc)


Gluco <- read.csv("A.Gluconeogenesis.csv", sep=",",header=F)

Glyco <- read.csv("B.Glyocgenolysis.csv", sep=",",header=F)

TCA <- read.csv("C.TCA.csv", sep=",",header=F)

Glyox <- read.csv("D.Glyoxylate.csv", sep=",",header=F)

Gluta <- read.csv("E.Glutamate.csv", sep=",",header=F)

Nucl <- read.csv("F.Nucleotide.csv", sep=",",header=F)

Amino <- read.csv("G.Amino.csv", sep=",",header=F)

Lipid <- read.csv("H.Lipid.csv", sep=",",header=F)


pdf(file="sampleCbModel_5000.pdf",width=30,height=60, pointsize=6);
lcex = 1.5;
par(mfrow=c(6,3),mar=c(15, 10, 25, 5));

z= seq(1,8);

for (i in 1:18) {

	f=cbind(Gluta[,i],TCA[,i],Glyox[,i],Gluco[,i],Glyco[,i],Nucl[,i],Amino[,i],Lipid[,i])

	boxplot(f,range=0,col='cyan',axes=FALSE)
	axis(1,at=z, labels=c("Glutamate formation","TCA cycle", "Glyoxylate cycle","Gluconeogenesis", "Glycogenolysis",  "Nucleotide", "Amino acid", "Lipid"),cex.axis=5,font=2,las = 2);
	axis(2,cex.axis=5,font=2,cex.lab=5);
}


dev.off();

