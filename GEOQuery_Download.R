###
### Script made Wed Nov  4 19:52:46 CST 2015
### R version 3.2.0 (2015-04-16) -- "Full of Ingredients"
###   Platform: x86_64-unknown-linux-gnu (64-bit)
###   GEOquery_2.34.0     Biobase_2.28.0      BiocGenerics_0.14.0
###   
#install GEOquery
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
###################################################################
# If you only have a GSE's list use the next chunk of code ########
###################################################################

GSMnamesFromVecGSEnames<- function(vecGSEnames){
  #This function recieves a vector of GSEnames
  # Build up  a list with equal length to the vector of GSEnames 
  GSEObjects<-list(length=length(vecGSEnames))
  #Fill thath list with GSEObjects
  for(i in 1:length(vecGSEnames)){
    GSEObjects[i]<-getGEO(vecGSEnames[i],GSEMatrix=FALSE)
    print("We have theGSE Object"); print(i)
  }
    #Call to the function to GEt a vector of GSM names of each GSE
    AllGSMnames<-NamesVectorGSMFromGSE(GSEObjects)
    #Call to a function to download the CEL files, from the previous GSMnames vector.
    DownCELFromThisGSMVec(AllGSMnames)
    #Return the GSMnames vector
  return(AllGSMnames)  
}
  # This function recieves a List of GSEObjects and return all GSMnames in a Vector 
NamesVectorGSMFromGSE <-function(ThisGSEObject){
 #Build up a vector of length zero
  vecpartialnames<-vector(length=0)
  #Fill the previous object with the GSM names.
    for (i in 1:length(ThisGSEObject)){
    vecpartialnames<-c(vecpartialnames,names(GSMList(ThisGSEObject[[i]])))
    }
  return(vecpartialnames)
}  

############################################################################  
#Once we have a GSMnames vector we can Download the .CEL files through:  ###
############################################################################
  
#Función calada para descargar los datos supplementarios (Incluyendo los CEL files) a partir de un vector de nombres GSM
#Para comenzar desde donde se quedo y evitar el problema es que las descargas se cortan  
DownCELFromThisGSMVec<-function(ThisGSMList){
    for(i in 1:5){
    sapply(ThisGSMList,function(x){ getGEOSuppFiles(x)})
    #ponle un if para comprobar exito en la descarga si no que comience
      if( length(ThisGSMList) == as.numeric(try(system("ls GSM*/*CEL* |  wc -l"))) ) {
      return()
      } 
    }
  }

#####################################################
###### More pseudo-random code in construction ######
#####################################################


  RobustAllGSMfromthisGSMList<-function(ThisGSMList){
    for(i in 1:20){
    sapply(ThisGSMList,function(x){ getGEOSuppFiles(x)})
    #ponle un if para comprobar exito en la descarga si no que comience
    
      if( length(ThisGSMList) == as.numeric(try(system("ls GSM*/*CEL* |  wc -l"))) ) {
      return()
      } 
      # Obten un archivo csv a partir de los nombres de los .CEL ya descargados
      (try(system("ls GSM*/*CEL* |  sed -e 's/\//_/g' already | sed 's/\./_/g' | sed 's/_/"\t"/g' | awk '{ print $2 }' | sed ':a;N;s/\n/,/g;ba' > CELAlreadyDown.csv"))
      #Convierte esa csv en un vector
      VecAlreadyDown<-as.vector(read.csv("CELAlreadyDown.csv"))
    }
  }


#Si no se descargo completamente
(names(GSMList(gse54002)))[!(names(GSMList(gse54002)) %in% vectc)]
sapply(as.list((names(GSMList(gse54002)))[!(names(GSMList(gse54002)) %in% vectc)]),function(x){ getGEOSuppFiles(x)})




DownCELFromThisGSMVec<-function(ThisGSMList){
    for(i in 1:5){
    sapply(ThisGSMList,function(x){ getGEOSuppFiles(x)})
    #ponle un if para comprobar exito en la descarga si no que comience
    try(system("ls GSM*/*CEL* |  wc -l > CELAlreadyDown.csv"))
    CelAD<-as.numeric(read.table("CELAlreadyDown.csv"))
      if( length(ThisGSMList) == CelAD ) {
      return()
      } 
    }
  }
  
  
  try(system("ls *CEL* |  wc -l > CELAlreadyDown.csv"))
  CelAD<-as.numeric(read.table("CELAlreadyDown.csv"))
    


