#Install and load the package
install.packages('ggfortify')
library(ggfortify)
# Define our group
myGSE= c(rep("GSE42568",121),rep("GSE50567",41),rep("GSE4002",433),rep("GSE10810",58),rep("GSE29431",66))
ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),rep("s",27),rep("t",31),rep("s",12),rep("t",54))
# We transpose the matrix 
RawExprs<-exprs(Data)
t.RawExprs<-t(RawExprs)
t.r.e<-cbind(t.RawExprs,myGSE)
t.r.e.L<-cbind(t.RawExprs,ListaSanosYEnf)
d.f.t.r.e <-as.data.frame(t.r.e)
d.f.t.r.e.L <-as.data.frame(t.r.e.L)
# Calculate the principal components and make the plot
pdf("PCA_Datos_Crudos_agrupados_por_GSE_SanosyEnfermos_719.pdf",width=7,height=5)
autoplot(prcomp(t.r.e), data=d.f.t.r.e, colour='myGSE', main="Principal components 719 raw data by GSE")
autoplot(prcomp(t.r.e.L), data=d.f.t.r.e.L, colour='ListaSanosYEnf', main="Principal components 719 raw data by health or sick")
dev.off()
2+1

ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),rep("s",27),rep("t",31),rep("s",12),rep("t",54))
# We transpose the matrix 
t.edata<-t(edata)
t.edata.ListaSanosYEnf<-cbind(t.edata,ListaSanosYEnf)
d.f.t.e <-as.data.frame(t.edata.ListaSanosYEnf)
# Calculate the principal components and make the plot
pdf("PCA_DatosNormalizados_only_fRMA719_Enfermos_y_Sanos.pdf",width=7,height=5)
autoplot(prcomp(t.edata), data=d.f.t.e, colour='ListaSanosYEnf', main="Principal components 719 only fRMA normalized")
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
GSM122nnnid<-names(GSM122nnnid)
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

ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),rep("s",27),rep("t",31),rep("s",12),rep("t",54))

GSM122nnnid<-read.csv("GSM122nnnid.txt")
GSM122nnnid<-names(GSM122nnnid)
