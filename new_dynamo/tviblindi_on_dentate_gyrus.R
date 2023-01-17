data<-readRDS("/data/results/dentategyrus/PCA_data_DentateGyrus.RDS")
group_id<-readRDS("/data/results/dentategyrus/cell_cycles_phase.RDS")
thematrrixumap<-readRDS("/data/results/dentategyrus/DentateGyrus_UMAP.RDS")
library(tviblindi)
tv1<-tviblindi(data=data, labels=group_id)
DimRed(tv1)
DimRed(tv1,method="umap")
DimRed(tv1,layout = thematrrixumap)
KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))


velocity_dentate <- readRDS("/data/results/dentategyrus/DentateGyrus_PCA_velocity.RDS")
#pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 30,MatrixOfVelocity = velocity_dentate)
new_pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 30,MatrixOfVelocity = velocity_dentate)
tv1$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv1$origin$calculatedPseudotimeNoRoot<-which.min(tv1$pseudotime$calculatedPseudotimeNoRoot$res)

tv1$pseudotime$calculatedPseudotimeMinRest$res<-as.numeric(new_pseudotime_no_root)
tv1$origin$calculatedPseudotimeMinRest<-which.min(tv1$pseudotime$calculatedPseudotimeMinRest$res)

Walks(tv1,N=1000,origin_name = "calculatedPseudotimeNoRoot")
Walks(tv1,N=1000,origin_name = "calculatedPseudotimeMinRest")
launch_shiny(tv1)
