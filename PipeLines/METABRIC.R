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

############################
#### BEGIN THE ANOTATION ###
############################

lumis_id_METABRIC<-rownames(Matrix.METABRIC)
write.table(lumis_id_METABRIC, file="lumis_id_METABRIC.txt", quote = F, row.names = F, col.names = F)

genesymbols<-getSYMBOL(as.character(lumis_id_METABRIC), "illuminaHumanv3.db")
write.table(genesymbols, file="illuminaHumanv3_geneSymbols.txt", quote = F, row.names = TRUE, col.names = F)




# Write genes 45 METABRIC  write.table(METABRIC_COLgenes45, file="METABRIC_COLgenes45.txt", quote = F, row.names = TRUE, col.names = T)



# My oun annotation of the chip 
# GenePROFILER to complete the poor annotation
# Perl for remove the spaces perl -pi.old -e 's{ }{-}g' Myannotation_illuminaHT12v3.csv
# my_annotation <- read.table(file="/home/rmejia/Documents/Doctorado/METABRIC/Myannotation_illuminaHT12v3NOspaces.txt")
# read_my_annotation = my_annotation 
# read_my_annotation.matrix<-as.matrix(read_my_annotation)



### Shell script pata quitar los primeros 48 mil y pico renglones  cat /home/rmejia/Documents/Doctorado/METABRIC/Myannotarion_illuminaHT12v3.txt  | awk ' NR >=48804   { print }' > /home/rmejia/Documents/Doctorado/METABRIC/Myannotation_illuminaHT12v3.txt
my_annotation <- as.matrix(read.table(file="/home/rmejia/Documents/Doctorado/METABRIC/Myannotation_illuminaHT12v3.txt", header=TRUE,header = FALSE, row.names=1,colClasses = "character"))



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

which(rownames(METABRIC_COL) %in% genes50)
which(rownames(METABRIC_COL) %in% genes50)

#faltan estos   "CDCA1" "CXXC5" "KNTC2" "MIA"   "ORC6L"  
#:O


cut -f 145-1141 coco.txt > cocoDisc.txt

 cat raul_metabric_pam50scores.txt  | awk ' NR <=2136   { print }' > pamscores_NONA.txt
which(rownames(RawExpMET) %in% colnames(Normexprs))
RawExpMET
groups=
WW.Matrix2<-as.numeric(RawExpMET[,1:5])
WW.Matrix2<-(matrix(as.numeric(RawExpMET[,1:5]),ncol=5))


myGSE= c(rep("N",143),rep("D",997),rep("V",995),"N")
ListaSanosYEnf<-c(rep("s",144),rep("t",1992))
LabelObject<-data.frame(myGSE,ListaSanosYEnf)
t.RawExprs<-WW.Matrix2
pdf("PCA_METABRICWrong.pdf",width=7,height=5)
autoplot(prcomp(t.RawExprs), data=LabelObject, colour='myGSE', main="PCA METABRIC Wrong")
dev.off()
1+2


