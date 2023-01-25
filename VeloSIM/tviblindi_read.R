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
simulated_data <- readRDS("/data/VeloSIM/datasets/2000_100_0.2.rds")

library(tviblindi)
#tv1<-tviblindi(data=scale(t(simulated_data$counts_s)),labels=as.factor(simulated_data$backbone))
tv1<-tviblindi(data=(simulated_data$counts_s_PCA),labels=as.factor(simulated_data$backbone))
DimRed(tv1)
DimRed(tv1,method="umap")

Set_origin(tv1,label = "5_6",origin_name = "5_6_hitting_time")

KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = FALSE,origin_name = "5_6_hitting_time")
Walks(tv1,N=1000,origin_name = "5_6_hitting_time")

pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 50, MatrixOfVelocity = simulated_data$velocity_100_pca)
tv1$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv1$origin$calculatedPseudotimeNoRoot<-which.min(tv1$pseudotime$calculatedPseudotimeNoRoot$res)
Walks(tv1,N=1000,origin_name = "calculatedPseudotimeNoRoot")

pseudotime_30 <- get_pseudotime_from_velocity(tv1, 30, IndexOfRootCell=tv1$origin$`5_6_hitting_time`,MatrixOfVelocity = t(simulated_data$velocity))
tv1$pseudotime$calculatedPseudotime30neighbours$res<-as.numeric(pseudotime_30)
tv1$origin$calculatedPseudotime30neighbours<-tv1$origin$`5_6_hitting_time`
Walks(tv1,N=1000,origin_name = "calculatedPseudotime30neighbours")


pseudotime_80 <- get_pseudotime_from_velocity(tv1, 80, IndexOfRootCell=tv1$origin$`5_6_hitting_time`,MatrixOfVelocity = t(simulated_data$velocity))
tv1$pseudotime$calculatedPseudotime80neighbours$res<-as.numeric(pseudotime_80)
tv1$origin$calculatedPseudotime80neighbours<-tv1$origin$`5_6_hitting_time`
Walks(tv1,N=1000, K=80, origin_name = "calculatedPseudotime80neighbours")


launch_shiny(tv1)
