

data<-readRDS("/data/results/Hematopoiesis/hematopoiesis_PCA.RDS")
group_id<-readRDS("/data/results/Hematopoiesis/hematopoiesis_cell_type.RDS")
thematrrixumap<-readRDS("/data/results/Hematopoiesis/hematopoiesis_UMAP.RDS")
library(tviblindi)
tv1<-tviblindi(data=data,labels=group_id)
DimRed(tv1)
DimRed(tv1,method="umap")
DimRed(tv1,layout = thematrrixumap)
# make layout col and rownames NULL
Set_origin(tv1,label = "HSC",origin_name = "HSC_hitting_time")
Set_origin(tv1,label = "HSC",origin_name = "HSC_hitting_distance")

KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = FALSE,origin_name = "HSC_hitting_time")
Walks(tv1,N=1000,origin_name = "HSC_hitting_time",K = 15)

Pseudotime(tv1,weighted = TRUE,origin_name = "HSC_hitting_distance") 
Walks(tv1,N=1000,origin_name = "HSC_hitting_distance") 

launch_shiny(tv1)

