.source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
#crear un objeto para asignar el org.Hs.egPATH2EG que mapea identificadores de genes entrez a los 
#identificadores utilizados por KEGG para las vías. 
kegg <- org.Hs.egPATH2EG
#asignar a un objeto el subconjunto de keys (x) mapeadas que fueron asignadas. 
mapped <- mappedkeys(kegg)
#Convertir el objeto en una lista. 
kegg2 <- as.list(kegg[mapped])

class(kegg)
class(org.Hs.egPATH2EG)
?mappedkeys
?org.Hs.egPATH2EG


kegg2
max(unlist(lapply(kegg2,length)))
hist(unlist(lapply(kegg2,length)))
?hist

biocLite("gage")
keeg_sets_uptodate<-kegg.gsets(species = "hsa", id.type = "kegg")
