

exprsFile="discovery_ExpressionMatrix.txt"




library("Biobase")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
minimalSet <- ExpressionSet(assayData=exprs)


exprsFileNorm="../Normals/normals_ExpressionMatrix.txt"
library("Biobase")
Normexprs <- as.matrix(read.table(exprsFileNorm, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
NormalminimalSet <- ExpressionSet(assayData=Normexprs)

proc.time()
exprsFileNorm="../Normals/normals_ExpressionMatrix.txt"
library("Biobase")
Normexprs <- as.matrix(read.table(exprsFileNorm, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
NormalminimalSet <- ExpressionSet(assayData=Normexprs)
