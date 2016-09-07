library(annotate)
library(hgu133plus2.db)

## Annotate a secas
#Ejemplo Hello<-select(hgu133plus2.db, "207356_at", c("SYMBOL","ENTREZID"))

path_to_probes<-c("/my/path/para/probes")
probesids<-read.table(path_to_probes)
SelectAffy<-select(hgu133plus2.db, as.character(probesids), c("SYMBOL","ENTREZID"))
write.table(SelectAffy,file="YourProbes_anoted")
