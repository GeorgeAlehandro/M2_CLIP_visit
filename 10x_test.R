data<-read.csv("/data/10X_analysis/X_pca.csv")
group_id<-tv0$labels[[1]]

library(tviblindi)
tv3<-tviblindi(data=as.matrix(data),labels=group_id)
DimRed(tv3)
DimRed(tv3,method="umap")

Set_origin(tv3,label = "A_B_blasts",origin_name = "A_B_blasts")


KNN(tv3)
Som(tv3,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv3) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv3,alpha2=1))

Pseudotime(tv3,weighted = TRUE,origin_name = "A_B_blasts") 
Walks(tv3,N=1000,origin_name = "A_B_blasts") 

tv3$vae$save("/data/10X_analysis/model.json", "/data/10X_analysis/model.h5")
##Adding un-rooted RNA velocity method
velocity <- read.csv("/data/10X_analysis/Velocity_pca.csv")

pseudotime_no_root <- get_pseudotime_from_velocity(tv3, 30, MatrixOfVelocity=as.matrix(velocity))
tv3$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv3$origin$calculatedPseudotimeNoRoot<-which.min(tv3$pseudotime$calculatedPseudotimeNoRoot$res)
Walks(tv3,N=1000,origin_name = "calculatedPseudotimeNoRoot")


singler_labels <- readRDS("10X_analysis/singleR_annotation.RDS")
tv3$labels<-singler_labels
launch_shiny(tv3)


