
listpathifier_to_matrix<-function(mylist){
  result<-mylist[[1]]
  for(i in 2:length(mylist)){
    result<-rbind(result,mylist[[i]])
  }
  rownames(result)<-names(mylist)
  return(result)
}

## Example of use
## If you have your list of PDS named "PDS$scores"
# matrixpathifier<-listpathifier_to_matrix(PDS$scores)
