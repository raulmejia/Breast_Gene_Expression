###
### Script made Wed Nov  4 19:52:46 CST 2015
### R version 3.2.0 (2015-04-16) -- "Full of Ingredients"
###   Platform: x86_64-unknown-linux-gnu (64-bit)
###   GEOquery_2.34.0     Biobase_2.28.0      BiocGenerics_0.14.0
###   
#install GEOquery

library(GEOquery)

# If you only have a GSE's list use the next chunk of code
gse54002 <- getGEO("GSE54002",GSEMatrix=FALSE)
sapply(as.list(names(GSMList(gse54002))),function(x){ getGEOSuppFiles(x)})





already<-list("GSM1305324","GSM1305364","GSM1305411","GSM1305517","GSM1305325","GSM1305365","GSM1305412","GSM1305518","GSM1305326","GSM1305366","GSM1305413","GSM1305519","GSM1305327","GSM1305367","GSM1305414","GSM1305520","GSM1305328","GSM1305368","GSM1305415","GSM1305521","GSM1305329","GSM1305369","GSM1305416","GSM1305522","GSM1305330","GSM1305370","GSM1305417","GSM1305523","GSM1305331","GSM1305371","GSM1305418","GSM1305524","GSM1305332","GSM1305372","GSM1305419","GSM1305539","GSM1305333","GSM1305373","GSM1305420","GSM1305540","GSM1305334","GSM1305374","GSM1305421","GSM1305541","GSM1305335","GSM1305375","GSM1305422","GSM1305542","GSM1305336","GSM1305376","GSM1305423","GSM1305543","GSM1305337","GSM1305377","GSM1305424","GSM1305544","GSM1305338","GSM1305378","GSM1305425","GSM1305545","GSM1305339","GSM1305379","GSM1305426","GSM1305546","GSM1305340","GSM1305380","GSM1305448","GSM1305547","GSM1305341","GSM1305381","GSM1305449","GSM1305548","GSM1305342","GSM1305382","GSM1305450","GSM1305549","GSM1305343","GSM1305383","GSM1305451","GSM1305550","GSM1305344","GSM1305384","GSM1305452","GSM1305551","GSM1305345","GSM1305385","GSM1305453","GSM1305552","GSM1305346","GSM1305386","GSM1305454","GSM1305577","GSM1305347","GSM1305387","GSM1305455","GSM1305578","GSM1305348","GSM1305388","GSM1305456","GSM1305579","GSM1305349","GSM1305389","GSM1305457","GSM1305580","GSM1305350","GSM1305390","GSM1305503","GSM1305581","GSM1305351","GSM1305391","GSM1305504","GSM1305582","GSM1305352","GSM1305392","GSM1305505","GSM1305583","GSM1305353","GSM1305393","GSM1305506","GSM1305584","GSM1305354","GSM1305394","GSM1305507","GSM1305663","GSM1305355","GSM1305402","GSM1305508","GSM1305719","GSM1305356","GSM1305403","GSM1305509","GSM1305723","GSM1305357","GSM1305404","GSM1305510","GSM1305728","GSM1305358","GSM1305405","GSM1305511","GSM1305729","GSM1305359","GSM1305406","GSM1305512","GSM1305730","GSM1305360","GSM1305407","GSM1305513","GSM1305731","GSM1305361","GSM1305408","GSM1305514","GSM1305732","GSM1305362","GSM1305409","GSM1305515","GSM1305733","GSM1305363","GSM1305410","GSM1305516","GSM1305734")
alreadyvector <- as.vector(already)
(names(GSMList(gse54002)))[!(names(GSMList(gse54002)) %in% already2)]
sapply(as.list((names(GSMList(gse54002)))[!(names(GSMList(gse54002)) %in% already2)]),function(x){ getGEOSuppFiles(x)})

#no se descargo completamente
(names(GSMList(gse54002)))[!(names(GSMList(gse54002)) %in% vectc)]
sapply(as.list((names(GSMList(gse54002)))[!(names(GSMList(gse54002)) %in% vectc)]),function(x){ getGEOSuppFiles(x)})


#Automatic show

DownObjectsGSE<- function(ListaGSEs){
GSEObjects<-list(length=length(ListaGSEs))
  for(i in 1:length(ListaGSEs)){
    GSEObjects[i]<-getGEO(ListaGSEs[i],GSEMatrix=FALSE)
    print("We have the first GSE Object")
    }
  return(GSEObjects)  
  }
  
  ALLGSMfromthisGSEObject<-function(ThisGSEObject){
    if()
    sapply(as.list(names(GSMList(ThisGSEObject))),function(x){ getGEOSuppFiles(x)})
    #ponle un if para comprobar exito en la descarga si no que comience
  }

  ALLGSMfromthisGSMList<-function(ThisGSMList){
    for(i in 1:20){
    sapply(ThisGSMList,function(x){ getGEOSuppFiles(x)})
    #ponle un if para comprobar exito en la descarga si no que comience
      if( length(ThisGSMList) == NumberAlreadyDown<-try(system("ls GSM* | grep *CEL* | wc -l")) ) {
      return()
      }
    }
  }

# If you have the exact list of GSM, you can use the next code:
# First we create a list object in R with the names of our GSM
MyGSMList <- list('','','')
# Next we download the suppl files (wich include the .CEL files).  The function will create a folder for 
# each GSM information.
sapply(MyGSMList,function(x){ getGEOSuppFiles(x)})

#We must make the list of GSM to download:

#Order the download

#mv GSM*/*.CEL* .

