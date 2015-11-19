#source("http://bioconductor.org/biocLite.R")
#biocLite("hgu133plus2.db")
#biocLite("hgu133plus2frmavecs")
#library(hgu133a.db) Instala las apropiadas

library(affy)
library(frma)
library(sva)
library(annotate)  
library(limma)
library(hgu133plus2.db)
library(hgu133plus2frmavecs)

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
mycolors= c(rep("purple",121),rep("yellow",41),rep("red",433),rep("green",58),rep("brown",66))

pdf("DatosCrudos719.pdf",width=7,height=5)
hist(Data, col=mycolors, main="Raw data distribution")
boxplot(Data,col=mycolors, main="Raw data distribution")

plot(100*pm.mm, type='h', main='Percent of MMs > PMs', ylab="%",xlab="Microarrays", ylim=c(0,50), col="red", lwd=5 )
grid(nx = NULL, ny = 6, col = "blue", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)

dev.off()

frmaData <- frma(Data, summarize="robust_weighted_average")
edata<-exprs(frmaData)

pdf("frmaNormalized.pdf",width=7,height=5)
plotDensity(edata, col=mycolors, main="frma normalization")
boxplot(edata,col=mycolors, main="Normaliced data distribution")
dev.off()

################# BATCH #########################

nombresGSM <- colnames(edata)

nombresGSM <- data.frame(lapply(nombresGSM, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

nombresGSM <- as.character(nombresGSM)

colnames(edata) <- nombresGSM

batch <- c((rep(0,719)))

# Enfermos

GSE1456<-read.table(/"Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE1456.txt", colClasses = "character")
GSE1561<-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE1561.txt", colClasses = "character")
GSE2603<-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE2603.txt", colClasses = "character")
GSE2990<-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE2990.txt", colClasses = "character")
GSE3494<-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE3494.txt", colClasses = "character")
GSE4922<-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE4922.txt", colClasses = "character")
GSE7390<-read.table("/Users/hachepunto/notron_media/mega_download/lists_of_names/all_GSE7390.txt", colClasses = "character")

GSE1456 <- GSE1456[[1]]
GSE1561 <- GSE1561[[1]]
GSE2603 <- GSE2603[[1]]
GSE2990 <- GSE2990[[1]]
GSE3494 <- GSE3494[[1]]
GSE4922 <- GSE4922[[1]]
GSE7390 <- GSE7390[[1]]

n1456 <- which(nombresGSM %in% GSE1456)
n1561 <- which(nombresGSM %in% GSE1561)
n2603 <- which(nombresGSM %in% GSE2603)
n2990 <- which(nombresGSM %in% GSE2990)
n3494 <- which(nombresGSM %in% GSE3494)
n4922 <- which(nombresGSM %in% GSE4922)
n7390 <- which(nombresGSM %in% GSE7390)

batch[n1456] = 1
batch[n1561] = 2
batch[n2603] = 3
batch[n2990] = 4
batch[n3494] = 5
batch[n4922] = 6
batch[n7390] = 7

tumor <- c(n1456, n1561, n2603, n2990, n3494, n4922, n7390)

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









