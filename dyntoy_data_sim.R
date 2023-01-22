devtools::install_github('dynverse/dyntoy')
devtools::install_github("dynverse/dyno")
library(dyno)
library(tidyverse)
n_events <- 5000
n_features <- 100
set.seed(12345)
d <- dyntoy::generate_dataset(
  id           = 'tviblindi_dyntoy_test',
  model        = 'connected',
  num_features = n_features,
  num_cells    = n_events,
  add_velocity = T
)

# velocity data stored at: d$expression_projected

data<-as.matrix(d[["expression"]])
group_id<-d[["progressions"]][["from"]]

library(tviblindi)
tv1_dyno_data<-tviblindi(data=data,labels=group_id)
DimRed(tv1_dyno_data)
DimRed(tv1_dyno_data,method="umap")
simplify_trajectory(tv1_dyno_data)
Set_origin(tv1_dyno_data,label = "M4",origin_name = "M4_hitting_time")
Set_origin(tv1_dyno_data,label = "M4",origin_name = "M4_hitting_distance")

KNN(tv1_dyno_data)
Som(tv1_dyno_data,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1_dyno_data) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1_dyno_data,alpha2=1))

Pseudotime(tv1_dyno_data,weighted = FALSE,origin_name = "M4_hitting_time")
Walks(tv1_dyno_data,N=1000,origin_name = "M4_hitting_time")

Pseudotime(tv1_dyno_data,weighted = TRUE,origin_name = "M4_hitting_distance") 
Walks(tv1_dyno_data,N=1000,origin_name = "M4_hitting_distance") 

launch_shiny(tv1_dyno_data)