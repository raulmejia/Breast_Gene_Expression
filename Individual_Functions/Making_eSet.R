

exprsFile=".txt"
library("Biobase")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
minimalSet <- ExpressionSet(assayData=exprs)

###############

exprsFileNorm="../Normals/normals_ExpressionMatrix.txt"
library("Biobase")
Normexprs <- as.matrix(read.table(exprsFileNorm, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
NormalminimalSet <- ExpressionSet(assayData=Normexprs)

proc.time()

exprsFileDisc="../DiscoverySet_Complete_4800/discovery_ExpressionMatrix.txt"
library("Biobase")
Discexprs<- as.matrix(read.table(exprsFileDisc, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
DiscminimalSet <- ExpressionSet(assayData=Discexprs)


exprsFileVal="../ValidationSet/validation_ExpressionMatrix.txt"
Valexprs<- as.matrix(read.table(exprsFileVal, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
ValminimalSet <- ExpressionSet(assayData=Valexprs)


Matrix.METABRIC<-cbind(Normexprs,Discexprs,Valexprs)
eSet.METABRIC<-ExpressionSet(assayData=Matrix.METABRIC)

###  To Do
# Read the data with more significant digits
