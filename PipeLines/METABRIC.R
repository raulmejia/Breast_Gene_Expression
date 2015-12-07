source("http://bioconductor.org/biocLite.R")
biocLite("illuminaHumanv3.db")

Initial.time=proc.time()

library(illuminaHumanv3.db)
library(limma)
library(annotate)
library(Biobase)

Initial.time=proc.time()

exprsFileNorm="../Normals/normals_ExpressionMatrix.txt"
Normexprs <- as.matrix(read.table(exprsFileNorm, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
NormalminimalSet <- ExpressionSet(assayData=Normexprs)


exprsFileDisc="../DiscoverySet_Complete_4800/discovery_ExpressionMatrix.txt"
Discexprs<- as.matrix(read.table(exprsFileDisc, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
DiscminimalSet <- ExpressionSet(assayData=Discexprs)


exprsFileVal="../ValidationSet/validation_ExpressionMatrix.txt"
Valexprs<- as.matrix(read.table(exprsFileVal, header=TRUE))
ValminimalSet <- ExpressionSet(assayData=Valexprs)

Matrix.METABRIC<-cbind(Normexprs,Discexprs,Valexprs)
eSet.METABRIC<-ExpressionSet(assayData=Matrix.METABRIC)

#Design and contingence matrix
design = matrix(rep(0,4272), nrow=2136)
colnames(design) = c('tumor','sano')
rownames(design) = colnames(Matrix.METABRIC)
design[145:2136,1]=1
design[1:144,2]=1
cont.matrix = makeContrasts('tumor - sano', levels=design)

fit = lmFit(Matrix.METABRIC, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

topTable(fit2, coef=1, adjust='fdr')

statistics <- topTable(fit2, coef=1, adjust='fdr', n = length(row.names(Matrix.METABRIC)), sort = "none")

summary(statistics[,6])

# Control plot for log Fold Change

pdf(file="MatrixMETABRIClogFC_statisitics.pdf")
boxplot(statistics[,1])
dev.off()

# Control plot for B statistic
pdf(file="MatrixMETABRIC_B_statisitics.pdf")
boxplot(statistics[,6])
dev.off()

lumis_id_METABRIC<-rownames(Matrix.METABRIC)
write.table(lumis_id_METABRIC, file="lumis_id_METABRIC.txt", quote = F, row.names = F, col.names = F)




genesymbols<-getSYMBOL(as.character(lumis_id_METABRIC), "illuminaHumanv3.db")
write.table(genesymbols, file="illuminaHumanv3_geneSymbols.txt", quote = F, row.names = TRUE, col.names = F)

# Read my oun annotation of the chip 
#my_annotation <- read.table(file="myannotation_plus2.txt", header = TRUE, row.names=1,colClasses = "character")

#precolaps_combat <- cbind(my_annotation, statistics[,6],Matrix.METABRIC)
#colnames(precolaps_combat)[2] <- c("b")
#write.table(precolaps_combat, file="METABRIC_matrix_precolaps.txt", quote = FALSE, sep = "\t", row.names = FALSE)


Final.time=proc.time() - Initial.time


any((as.vector(attributes(genesymbols)))[[1]] == rownames(Matrix.METABRIC))

M.METABRIC.biocAnot<-Matrix.METABRIC
rownames(M.METABRIC.biocAnot)<-genesymbols

precolaps_METABRIC <- cbind(genesymbols, statistics[,6],Matrix.METABRIC)
colnames(precolaps_METABRIC)[2] <- c("b")
write.table(precolaps_METABRIC, file="METABRIC_matrix_precolaps.txt", quote = FALSE, sep = "\t", row.names = FALSE)



### Colapsaitor..

## Shell script pata quitar el NA del renglon  cat METABRIC_matrix_precolaps_colapsed.txt  | awk 'NR <= 2 || NR >=4   { print }' > METABRIC_colapsed_NO_NA.txt
exprsFileDisc="METABRIC_colapsed_NO_NA.txt"
METABRIC_Colapsed<- as.matrix(read.table(exprsFileDisc, header=TRUE, sep="\t", row.names=1, as.is=TRUE))

##Eliminar NA adentro de la matrix (convertirlos acero)

Whichones<-is.finite(METABRIC_Colapsed_NOinternalNA)

which(colSums(Whichones) %in% "19306")
which(rowSums(Whichones) %in% "2134")
which(rowSums(Whichones) %in% "2135")
which(rowSums(Whichones) %in% "2136")
METABRIC_Colapsed[1509,1462]
## Meterles ceros a los NA
# Ejemplo m <- matrix(c(1,2,3,NA,4,5), 3); m[!is.finite(m)] <- 0
METABRIC_Colapsed_NOinternalNA[!is.finite(METABRIC_Colapsed_NOinternalNA)] <- 0


# NA errased METColNAerased[, !colSums(is.na(METColNAerased))]
na.omit()

biocLite("genefu")
