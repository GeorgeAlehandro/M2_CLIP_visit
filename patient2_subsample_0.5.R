data<-readRDS("/data/patient2_0.5_subsample/X_pca.RDS")
group_id<-readRDS("/data/patient2_0.5_subsample/group_id.RDS")
cell_barcodes <- read.csv("/data/patient2_0.5_subsample/dynamo_cell_barcodes.csv")
singler_cell_types <- read.csv("/data/patient2_0.5_subsample/singleR_cell_types.tsv", sep="\t")
gating_labels <- unlist(tv0$labels)[match(cell_barcodes[,1],fcs@exprs[,6])]
ordered_singler_cell_types <- singler_cell_types[match(cell_barcodes[,1],singler_cell_types$Cell_CB),]
tv1$labels$'SingleR' <- as.factor(ordered_singler_cell_types[,2])
tv1$labels$default = as.factor(gating_labels)
library(tviblindi)
tv1<-tviblindi(data=as.matrix(data),labels=gating_labels)
DimRed(tv1)
DimRed(tv1,method="umap")

Set_origin(tv1,label = as.integer(1955),origin_name = "Origin1955")
Pseudotime(tv1,weighted = TRUE,origin_name = "Origin1955") 
Walks(tv1,N=1000,origin_name = "Origin1955") 

Set_origin(tv1,label = "A_B-blasts",origin_name = "A_B-blasts")

KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = TRUE,origin_name = "A_B-blasts") 
Walks(tv1,N=1000,origin_name = "A_B-blasts") 

tv1$vae$save("/data/patient2_0.5_subsample/patient2.json", "/data/patient2_0.5_subsample/patient2.h5")
##Adding un-rooted RNA velocity method
velocity <- readRDS("/data/patient2_0.5_subsample/PCA_velocity.RDS")
pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 30, MatrixOfVelocity=velocity)
tv1$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv1$origin$calculatedPseudotimeNoRoot<-which.min(tv1$pseudotime$calculatedPseudotimeNoRoot$res)
Walks(tv1,N=1000,origin_name = "calculatedPseudotimeNoRoot")

pseudotime_50 <- get_pseudotime_from_velocity(tv1, 50, MatrixOfVelocity=velocity, tv1$origin$`A_B-blasts`)
tv1$pseudotime$calculatedPseudotime50neighbours$res<-as.numeric(pseudotime_50)
tv1$origin$calculatedPseudotime50neighbours<-tv1$origin$`A_B-blasts`
Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours")

pseudotime_50_1955 <- get_pseudotime_from_velocity(tv1, 50, MatrixOfVelocity=velocity, as.integer(1955))
tv1$pseudotime$calculatedPseudotime50neighbours1955$res<-as.numeric(pseudotime_50_1955)
tv1$origin$calculatedPseudotime50neighbours1955<-as.integer(1955)
Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours1955")

launch_shiny(tv1)


