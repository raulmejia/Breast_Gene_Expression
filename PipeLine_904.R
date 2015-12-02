#source("http://bioconductor.org/biocLite.R")
#biocLite("hgu133plus2.db")
#biocLite("hgu133plus2frmavecs")
#library(hgu133a.db) Instala las apropiadas
biocLite("hgu133plus2.db")
biocLite("hgu133plus2frmavecs")
biocLite("limma")
biocLite("sva")
biocLite("frma")

# install.packages("gsubfn)"

#Nota: Haz un pca para cada etapa de los datos, tambien puede ser un herarchical clustering. para ver que tanto
#hemos combatido el batch effect

library(affy)
library(frma)
library(sva)
library(annotate)  
library(limma)
library(hgu133plus2.db)
library(hgu133plus2frmavecs)
library(ggfortify)
#library(gsubfn)

Data<-ReadAffy()
N=length(Data@phenoData@data$sample)
pm.mm=0
for (i in 1:N) {pm.mm[i] = mean(mm(Data[,i])>pm(Data[,i]))}
#GSE42568 (GSM1045*)=121
#GSE50567 (GSM1223*)=41
#GSE54002 (GSM1305*)=433
#GSE10810 (GSM272*)=58
###########GSE29044 (GSM719*)=79
#GSE29431 (GSM728)=66
#mycolors= c(rep(1,121),rep(2,41),rep(3,433),rep(4,58),rep(5,79),rep(6,66))
mycolors= c(rep("purple",121),rep("yellow",41),rep("red",433),rep("pink",185),rep("green",58),rep("brown",66))

pdf("DatosCrudos.pdf",width=7,height=5)
hist(Data, col=mycolors, main="Raw data distribution")
boxplot(Data,col=mycolors, main="Raw data distribution")

plot(100*pm.mm, type='h', main='Percent of MMs > PMs', ylab="%",xlab="Microarrays", ylim=c(0,50), col="red", lwd=5 )
grid(nx = NULL, ny = 6, col = "blue", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)

dev.off()

### The PCA
# We extract the expression matrix and transpose it to obtain the principal components 
RawExprs<-exprs(Data)
Inicial.time<-proc.time()
t.RawExprs<-t(RawExprs)
#make the Labels
myGSE= c(rep("GSE42568",121),rep("GSE50567",41),rep("GSE4002",433),rep("GSE10780",185),rep("GSE10810",58),rep("GSE29431",66))
ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),"s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","t","s","s","s","s","s","s","s","s","s","s","s","t","t","t","s","s","s","s","s","t","s","s","s","s","t","t","s","s","s","s","t","t","t","t","t","t","s","s","s","s","s","s","s","s","s","s","t","s","s","t","s","s","s","s","s","t","s","t","s","t","s","s","t","t","s","s","s","t","s","s","s","t","s","s","s","s","s","t","s","s","s","s","s","s","s","s","s","s","t","t","s","t","t","t","s","s","s","s","t","t","t","s","s","s","s","s","s","s","t","s","s","s","s","s","s","s","s","s","t","s","s","s","s","s","s","s","t","s","s","s","s","s","s","t","s","s","t","s","s","s","s","s","s","s","s","s","t","t","s","s","t","s","s","s","s","t","t","s","s","s","s","s","t","s",rep("s",27),rep("t",31),rep("s",12),rep("t",54))
PCARawLabels<-data.frame(myGSE,ListaSanosYEnf)

pdf("PCA_Datos_Crudos_GSE_colors_and_HealthvsSick_shape_719.pdf",width=7,height=5)
autoplot(prcomp(t.RawExprs), data=PCARawLabels, colour='myGSE',shape='ListaSanosYEnf', main="PCA 7904 raw data GSE color and Health vs Sick shape")
autoplot(prcomp(t.RawExprs), data=PCARawLabels, colour='ListaSanosYEnf',shape='myGSE', main="PCA 904 raw data GSE shape and Health vs Sick colour")
dev.off()
final.time<-proc.time() - Inicial.time
1+2
### End PCA





frmaData <- frma(Data, summarize="robust_weighted_average")
edata<-exprs(frmaData)

pdf("frmaNormalized.pdf",width=7,height=5)
plotDensity(edata, col=mycolors, main="frma normalization")
boxplot(edata,col=mycolors, main="Normalized data distribution")
dev.off()

Inicial.time<-proc.time()
t.RawExprs<-t(edata)
#make the Labels
myGSE= c(rep("GSE42568",121),rep("GSE50567",41),rep("GSE4002",433),rep("GSE10780",185),rep("GSE10810",58),rep("GSE29431",66))
ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),"s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","s","t","s","s","s","s","s","s","s","s","s","s","s","t","t","t","s","s","s","s","s","t","s","s","s","s","t","t","s","s","s","s","t","t","t","t","t","t","s","s","s","s","s","s","s","s","s","s","t","s","s","t","s","s","s","s","s","t","s","t","s","t","s","s","t","t","s","s","s","t","s","s","s","t","s","s","s","s","s","t","s","s","s","s","s","s","s","s","s","s","t","t","s","t","t","t","s","s","s","s","t","t","t","s","s","s","s","s","s","s","t","s","s","s","s","s","s","s","s","s","t","s","s","s","s","s","s","s","t","s","s","s","s","s","s","t","s","s","t","s","s","s","s","s","s","s","s","s","t","t","s","s","t","s","s","s","s","t","t","s","s","s","s","s","t","s",rep("s",27),rep("t",31),rep("s",12),rep("t",54))
PCARawLabels<-data.frame(myGSE,ListaSanosYEnf)

pdf("PCA_Datos_frmaOnlyGSE_colors_and_HealthvsSick_shape_719.pdf",width=7,height=5)
autoplot(prcomp(t.RawExprs), data=PCARawLabels, colour='myGSE',shape='ListaSanosYEnf', main="PCA 904 raw data GSE color and Health vs Sick shape")
autoplot(prcomp(t.RawExprs), data=PCARawLabels, colour='ListaSanosYEnf',shape='myGSE', main="PCA 904 raw data GSE shape and Health vs Sick colour")
dev.off()
final.time<-proc.time() - Inicial.time
1+2

################# BATCH #########################

nombresGSM <- colnames(edata)
nombresGSM <- data.frame(lapply(nombresGSM, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

#nombresGSM <- as.character(nombresGSM)

#colnames(edata) <- nombresGSM

batch <- c((rep(0,719)))

# Enfermos
#The next files must be in .csv format preferently
# Experimento con muchos controles GSE10780<-read.csv("/home/rmejia/Documents/Doctorado/Determinar_my_Pipeline_RM_Doctorado/GSE10780id.txt")

GSE54002<-read.csv("/home/rmejia/Documents/Doctorado/Determinar_my_Pipeline_RM_Doctorado/GSE54002id.txt")
GSE50567<-read.csv("/home/rmejia/Documents/Doctorado/Determinar_my_Pipeline_RM_Doctorado/GSE50567id.txt")
GSE42568<-read.csv("/home/rmejia/Documents/Doctorado/Determinar_my_Pipeline_RM_Doctorado/GSE42568id.txt")
GSE29431<-read.csv("/home/rmejia/Documents/Doctorado/Determinar_my_Pipeline_RM_Doctorado/GSE29431id.txt")
GSE10810<-read.csv("/home/rmejia/Documents/Doctorado/Determinar_my_Pipeline_RM_Doctorado/GSE10810id.txt")

GSE54002<-names(GSE54002)
GSE50567<-names(GSE50567)
GSE42568<-names(GSE42568)
GSE29431<-names(GSE29431)
GSE10810<-names(GSE10810)

n1305 <- which(nombresGSM %in% GSE54002)
n1223 <- which(nombresGSM %in% GSE50567)
n1045 <- which(nombresGSM %in% GSE42568)
n7286 <- which(nombresGSM %in% GSE29431)
n2729 <- which(nombresGSM %in% GSE10810)

batch[n1305] = 1
batch[n1223] = 2
batch[n1045] = 3
batch[n7286] = 4
batch[n2729] = 5

tumor <- c(n1305, n1223, n1045, n7286, n2729)

# Sanos

GSE15852 <-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE15852.txt", colClasses = "character")
GSE6883 <-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE6883.txt", colClasses = "character")
GSE9574 <-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE9574.txt", colClasses = "character")

GSE15852 <- GSE15852[[1]]
GSE6883 <- GSE6883 [[1]]
GSE9574 <- GSE9574 [[1]]

n15852 <- which(nombresGSM %in% GSE15852)
n6883  <- which(nombresGSM %in% GSE6883)
n9574  <- which(nombresGSM %in% GSE9574)

batch[n15852] = 8
batch[n6883] = 9
batch[n9574] = 10

healthy <- c(n15852, n6883, n9574)

sample <- c(1:880)
outcome <-c(rep("a",880))
outcome[tumor] <- "tumor"
outcome[healthy] <- "healthy"
cancer <- c(rep("a",880))
cancer[healthy] <- "normal"
pheno <- data.frame(sample,outcome, batch, cancer)
rownames(pheno) <- nombresGSM


group <- c(1:880)
group[tumor] = 1
group[healthy] = 2
group <- as.factor(group)

modcombat = model.matrix(~ group)

options("contrasts")
modcombat = model.matrix(~ outcome, pheno, contrasts = list(outcome = "contr.sum"))

head(modcombat)

modcombat[tumor,"outcome1"]

modcombat = model.matrix(~1, data=pheno)

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat)

write.csv(combat_edata, file="expression_matrix_affy.csv")

pdf("ComBat_Normalized.pdf",width=7,height=5)

mycolors = rep(c("blue","red","green", "magenta"), each = 2)
plotDensity(combat_edata, col=mycolors, main="ComBat normalization")
boxplot(combat_edata,col=mycolors, main="Normalized data distribution")
dev.off()


# Design and contingence matrix
design = matrix(rep(0,1760), nrow=880)
colnames(design) = c('tumor','sano')
rownames(design) = colnames(combat_edata)
design[n1456,1]=1
design[n1561,1]=1
design[n2603,1]=1
design[n2990,1]=1
design[n3494,1]=1
design[n4922,1]=1
design[n7390,1]=1
design[n15852,2]=1
design[n6883,2]=1
design[n9574,2]=1
cont.matrix = makeContrasts('tumor - sano', levels=design)

fit = lmFit(combat_edata, design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

topTable(fit2, coef=1, adjust='fdr')

statistics <- topTable(fit2, coef=1, adjust='fdr', n = length(row.names(combat_edata)), sort = "none")

head(statistics[,6])

write.csv(statistics, file = "statistics.csv", row.names=T, col.names=T)

affys<-rownames(combat_edata)

write.table(affys, file="affysIDs.txt", quote = F, row.names = F, col.names = F)

genesymbols<-getSYMBOL(as.character(affys), "hgu133a.db")
eset<-cbind(genesymbols, combat_edata)
eset<-ifelse(is.na(eset), as.vector(rownames(eset)), as.vector(eset))

write.csv(eset, file="precolaps.csv")
#####PCA DE Combat ###


t.matrix<-t(combat_edata)
#make the Labels
myGSE= c(rep("GSE42568",121),rep("GSE50567",41),rep("GSE4002",433),rep("GSE10810",58),rep("GSE29431",66))
ListaSanosYEnf<-c(rep("s",17),rep("t",104),rep("t",35),rep("s",6),rep("t",300),rep("s",16),rep("t",117),rep("s",27),rep("t",31),rep("s",12),rep("t",54))
PCARawLabels<-data.frame(myGSE,ListaSanosYEnf)

pdf("PCA_combat_GSE_colors_and_HealthvsSick_shape_719.pdf",width=7,height=5)
autoplot(prcomp(t.matrix), data=PCARawLabels, colour='myGSE',shape='ListaSanosYEnf', main="PCA 719 combat GSE color and Health vs Sick shape")
autoplot(prcomp(t.matrix), data=PCARawLabels, colour='ListaSanosYEnf',shape='myGSE', main="PCA 719 combat GSE shape and Health vs Sick colour")
dev.off()
1+2

### NEW Hugo

/home/hachepunto/rauldb/rauldb.R


### My Batch ~
mod=model.matrix(~as.factor(cancer) + as.factor(batch), data=pheno)
fit=lm.fit(mod,t(edata))
hist
