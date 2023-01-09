##first 25 principal components of dataset generated as follows

# devtools::install_github('dynverse/dyntoy')
# 
# n_events <- 10000
# n_features <- 3896 
# set.seed(12345)
# d <- dyntoy::generate_dataset(
#   id           = 'tviblindi_dyntoy_test',
#   model        = 'connected',
#   num_features = n_features,
#   num_cells    = n_events
# )


data<-read.csv("/data/X_pca.csv")
group_id<-read.csv("/data/cell_cycle_phase.csv")[2]

library(tviblindi)
tv1<-tviblindi(data=data,labels=group_id)
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
