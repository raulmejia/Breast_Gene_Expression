args <- commandArgs(trailingOnly = TRUE)
Path_to_your_matrix<-args[1]
#Path_to_your_matrix <-c("/castle/rmejia/Chucho_Biomarkers/HClustering/Data/OnlyProgression.txt")
#labels_path<-c("/castle/rmejia/Chucho_Biomarkers/HClustering/Data/CoLoRs_Labels_Controls_and_Normal_separated_TCGA.txt")
labels_path<-args[2]
#Results_path<-c("/castle/rmejia/Chucho_Biomarkers/HClustering/Results/")
Results_path<-args[3]
#names_to_save<-c("OnlyProg")
names_to_save<-args[4]
####################################
###### Loading and/or installing needig libraries
####################################
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("dendextend")) {
  install.packages("dendextend", dependencies = TRUE)
  library(dendextend)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

####################################
###### Reading the data
####################################

df <- read.table(Path_to_your_matrix)
labels <- read.table(labels_path)

all(colnames(df) == rownames(labels))
df<-(as.numeric(as.matrix(df)))

##################################
#### clustering and Visualizing
################################### 
#grp <- c(rep("green",112),rep("red",808))
dend <- df %>%  dist(method= "euclidean") %>% 
  hclust(method = "ward.D2") %>% as.dendrogram %>%
  set("branches_k_color", k=4) %>% set("branches_lwd", 1.2) %>%
  #set("labels_colors",grp) %>%
  set("labels_cex", c(.9,.9)) %>% 
  set("leaves_pch", 19) %>% set("leaves_col", labels[,1])
# plot the dend in usual "base" plotting engine:
plot(dend)
#the_bars <- cbind(grp, k_4d)
#dend %>% set("labels", "") %>% plot
#colored_bars(colors = the_bars, dend = dend)

ggd1 <- as.ggdend(dend)
ggplot(ggd1, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")

###############################
###### Visualizing
#################################

png(paste0(Results_path,names_to_save,c("_Rectangle"),c(".png")), # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # font size
plot(dend)
dev.off()

png(paste0(Results_path,names_to_save,c("_Circle"),c(".png")), # Name of png file
    width = 6 * 500,      # Easier scaling 6*500 = 3000 pixels
    height = 6 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # font size
ggplot(ggd1, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")
dev.off()


############################
### Getting the subgroups


###### 

for(w in 2:10){
  ko_as_d <- cutree(dend,k = w, order_clusters_as_data = TRUE) 
  df_k<-data.frame(ko_as_d,labels[,1])
  rownames(df_k)<-rownames(labels)
  write.table(df_k,file = paste0(Results_path,names_to_save,c("dend_"),w,c("_groups.txt")),sep="\t",quote = FALSE)
}

