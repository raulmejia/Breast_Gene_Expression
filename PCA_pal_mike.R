#Install and load the package
install.packages('ggfortify')
library(ggfortify)
# Define our group
myGSE= c(rep("GSE42568",121),rep("GSE50567",41),rep("GSE4002",433),rep("GSE10810",58),rep("GSE29431",66))
# We transpose the matrix 
t.edata<-t(edata)
t.edata.myGSE<-cbind(t.edata,myGSE)
d.f.t.e <-as.data.frame(t.edata.myGSE)
# Calculate the principal components and make the plot
pdf("PCA_DatosNormalizados_only_fRMA719.pdf",width=7,height=5)
autoplot(prcomp(t.edata), data=d.f.t.e, colour='myGSE', main="Principal components 719 only fRMA normalized")
dev.off()
