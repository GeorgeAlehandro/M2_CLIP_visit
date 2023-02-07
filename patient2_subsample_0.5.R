data<-readRDS("/data/patient2_0.5_subsample/X_pca.RDS")
group_id<-readRDS("/data/patient2_0.5_subsample/group_id.RDS")

library(tviblindi)
tv1<-tviblindi(data=as.matrix(data),labels=group_id)
DimRed(tv1)
DimRed(tv1,method="umap")

Set_origin(tv1,label = "M4",origin_name = "M4_hitting_time")
Set_origin(tv1,label = "M4",origin_name = "M4_hitting_distance")

KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = FALSE,origin_name = "M4_hitting_time")
Walks(tv1,N=1000,origin_name = "M4_hitting_time")

Pseudotime(tv1,weighted = TRUE,origin_name = "M4_hitting_distance") 
Walks(tv1,N=1000,origin_name = "M4_hitting_distance") 

launch_shiny(tv1)
