library(affy)
library(frma)
library(sva)
library(annotate)  
library(limma)
library(hgu133plus2.db)
library(hgu133plus2frmavecs)
library(ggfortify)


# Colors and groups for the graphics
`mycolors= c(rep("purple",121),rep("yellow",41),rep("red",433),rep("pink",185),rep("green",58),rep("brown",66))`
`myGSE= c(rep("GSE42568",121),rep("GSE50567",41),rep("GSE4002",433),rep("GSE10780",185),rep("GSE10810",58),rep("GSE29431",66))`
`ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),rep("s",27),rep("t",31),rep("s",12),rep("t",54))`

LabelObject<-data.frame(myGSE,ListaSanosYEnf)

### The PCA
# We extract the expression matrix and transpose it to obtain the principal components 

time.transpose_since<-proc.time()

myGSE= c(rep("Norm",144),rep("Disc",997),rep("Val",995))
ListaSanosYEnf<-c(rep("s",144),rep("t",1992))
LabelObject<-data.frame(myGSE,ListaSanosYEnf)
t.RawExprs<-t(RawExprs)
pdf("PCA_Datos_Crudos_GSE_colors_and_HealthvsSick_shape_719.pdf",width=7,height=5)
autoplot(prcomp(t.RawExprs), data=LabelObject, colour='myGSE',shape='ListaSanosYEnf', main="PCA 7904 raw data GSE color and Health vs Sick shape")
time.trasposeToPCA<-proc.time() - time.transpose_since
1+2

autoplot(prcomp(t.RawExprs), data=PCARawLabels, colour='ListaSanosYEnf',shape='myGSE', main="PCA 904 raw data GSE shape and Health vs Sick colour")
dev.off()
time.trasposeToPCA<-proc.time() - time.transpose_since
1+2
