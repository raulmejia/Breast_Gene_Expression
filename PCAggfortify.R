install.packages('ggfortify)
library(ggfortify)
my<-iris[c(1,2,3,4)]
autoplot(prcomp(my))

# We transpose the matrix 
t.edata<-t(edata)
t.edata.mycolors<-cbind(t.edata,mycolors)
d.f.t.e <-as.data.frame(t.edata.mycolors)
# Calculate the principal components and make the plot
autoplot(prcomp(t.edata), data=d.f.t.e, colour='mycolors')


autoplot(prcomp())
autoplot(prcomp(my), data=iris, colour='Species')


#GSE42568 (GSM1045*)=121
#GSE50567 (GSM1223*)=41
#GSE54002 (GSM1305*)=433
#GSE10810 (GSM272*)=58
###########GSE29044 (GSM719*)=79
#GSE29431 (GSM728)=66
