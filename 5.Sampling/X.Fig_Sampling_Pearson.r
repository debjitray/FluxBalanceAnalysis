library(Biobase)
library (Hmisc)


All <- read.csv("Pearson_5000.csv", sep=",",header=F)



pdf(file="sampleCbModel_Pears_5000.pdf",width=30,height=20, pointsize=6);

par(mfrow=c(1,1),mar=c(10, 10, 5, 2));


f=cbind(All[,1],All[,2],All[,3],All[,4],All[,5],All[,6],All[,7],All[,8],All[,9],All[,10],All[,11],All[,12],All[,13],All[,14],All[,15],All[,16],All[,17],All[,18]);

boxplot(f,range=0,col='cyan',axes=FALSE,ylim=c(-1,1))
axis(1,at=seq(1,18,1),labels=c(0.033,0.5,1,1.5,2,3,6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9.25,11),padj=1,font=2,cex.axis = 5);
axis(2,at=seq(-1,1,0.2),cex.axis=5,font=2,cex.lab=5);



dev.off();

