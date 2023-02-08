data<-readRDS("/data/patient2_0.5_subsample/X_pca.RDS")
group_id<-readRDS("/data/patient2_0.5_subsample/group_id.RDS")

library(tviblindi)
tv1<-tviblindi(data=as.matrix(data),labels=group_id)
DimRed(tv1)
DimRed(tv1,method="umap")

Set_origin(tv1,label = "B",origin_name = "M4_hitting_time")
Set_origin(tv1,label = "B",origin_name = "M4_hitting_distance")

KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = FALSE,origin_name = "M4_hitting_time")
Walks(tv1,N=1000,origin_name = "M4_hitting_time")

Pseudotime(tv1,weighted = TRUE,origin_name = "M4_hitting_distance") 
Walks(tv1,N=1000,origin_name = "M4_hitting_distance") 

tv1$vae$save("/data/patient2_0.5_subsample/patient2.json", "/data/patient2_0.5_subsample/patient2.h5")
##Adding un-rooted RNA velocity method
velocity <- readRDS("/data/patient2_0.5_subsample/PCA_velocity.RDS")
pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 30, MatrixOfVelocity=velocity)
tv1$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv1$origin$calculatedPseudotimeNoRoot<-which.min(tv1$pseudotime$calculatedPseudotimeNoRoot$res)
Walks(tv1,N=1000,origin_name = "calculatedPseudotimeNoRoot")

launch_shiny(tv1)
