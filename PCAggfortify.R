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


#GSE42568 (GSM1045*)=121
#GSE50567 (GSM1223*)=41
#GSE54002 (GSM1305*)=433
#GSE10810 (GSM272*)=58
###########GSE29044 (GSM719*)=79
#GSE29431 (GSM728)=66
MySanosVsSick<-c()


GSM104nnid<-read.csv("GSM104nnnid.txt")
GSM104nnid<-names(GSM104nnid)
GSM104nnsanos<-read.csv("GSM104nnnsanos.txt")
GSM104nnsanos<-names(GSM104nnsanos)

GSM122nnnid<-read.csv("GSM122nnnid.txt")
GSM122nnnid<-names(GSM122nnid)
GSM122nnnsanos<-read.csv("GSM122nnnsanos.txt")
GSM122nnnsanos<-names(GSM122nnnsanos)

GSM130nnnid<-read.csv("GSM130nnnid.txt")
GSM130nnnid<-names(GSM130nnnid)
GSM130nnnsanos<-read.csv("GSM130nnnsanos.txt")
GSM130nnnsanos<-names(GSM130nnnsanos)

GSM272o3nnnid<-read.csv("GSM272o3nnnid.txt")
GSM272o3nnnid<-names(GSM272o3nnnid)
GSM272o3nnnsanos<-read.csv("GSM272o3nnnsanos.txt")
GSM272o3nnnsanos<-names(GSM272o3nnnsanos)

GSM279nnnid<-read.csv("GSM279nnnid.txt")
GSM279nnnid<-names(GSM279nnnid)
GSM279nnnsanos<-read.csv("GSM279nnnsanos.txt")
GSM279nnnsanos<-names(GSM279nnnsanos)

GSM728nnnid<-read.csv("GSM728nnnid.txt")
GSM728nnnid<-names( GSM728nnnid)
GSM728nnnsanos<-read.csv("GSM728nnnsanos.txt")
GSM728nnnsanos<-names(GSM728nnnsanos)
