# Design and contingence matrix
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

genesymbols<-getSYMBOL(as.character(lumis_id_METABRIC), "hgu133plus2.db")
write.table(genesymbols, file="hgu133plus2DB.txt", quote = F, row.names = TRUE, col.names = F)

# Read my oun annotation of the chip 
my_annotation <- read.table(file="myannotation_plus2.txt", header = TRUE, row.names=1,colClasses = "character")

precolaps_combat <- cbind(my_annotation, statistics[,6],combat_edata)
colnames(precolaps_combat)[2] <- c("b")
write.table(precolaps_combat, file="rauldb_matrix_precolaps.txt", quote = FALSE, sep = "\t", row.names = FALSE)


source("http://biocondutor.org/biocLite.R")
biocLite("illuminaHumanv3.db")
library(illuminaHumanv3.db)










