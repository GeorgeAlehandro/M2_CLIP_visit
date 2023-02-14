## DONT FORGET MODIF_GATINGSET_TVi.R
fcstable("/data/fcs_ADT/filtered/filtered_adt.fcs",table_file = "/data/fcs_ADT/filtered/channels.txt")
tv0<-tviblindi_from_flowjo1(wsp = "/data/fcs_ADT/filtered/swALL_10x_gates_filtered.wsp",
                            fcsn = "/data/fcs_ADT/filtered/filtered_adt.fcs",eventsufix = "swALL$",
                            fcstable = "/data/fcs_ADT/filtered/channels.txt",origin = NULL,eventprefix = "[A-C]_.*")

###########
dynamo_filtered <- read.csv("/data/10X_analysis/gated_barcodes_dynamo_available.csv")$X0
seurat_expression_matrix <- read.csv("/data/10X_analysis/gene_expression_matrix.csv")
dyna1<-data.frame(barcode=dynamo_filtered,dynamo_order=1:length(dynamo_filtered))
seur1<-data.frame(barcode=gsub("\\.1$","",colnames(seurat_expression_matrix)[-1]),seurat_order=1:ncol(seurat_expression_matrix[,-1]))
both<-merge(seur1,dyna1,all.y = TRUE,all.x=FALSE)
#########
write.csv(tv0$events_sel,"/data/fcs_ADT/filtered/events_selected.csv")

indexes_to_keep <- readRDS(file = "10X_analysis/indexes_to_keep_after_dynamo_seurat.rds")
tv0$labels[[1]][indexes_to_keep]

data1<-readRDS("/data/10X_analysis/pca_embed_after_seurat.RDS")
group_id<-tv0$labels[[1]]
library(tviblindi)
tv1<-tviblindi(data=as.matrix(data1),labels=group_id)
DimRed(tv1)
DimRed(tv1,method="umap")

Set_origin(tv1,label = "A_B_blasts",origin_name = "A_B_blasts")


KNN(tv1)
Som(tv1,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv1) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = TRUE,origin_name = "A_B_blasts") 
Walks(tv1,N=1000,origin_name = "A_B_blasts") 

#tv1$vae$save("/data/patient2_0.5_subsample/patient2.json", "/data/patient2_0.5_subsample/patient2.h5")
##Adding un-rooted RNA velocity method
velocity <- readRDS("/data/10X_analysis/10x_PCA_velocity.RDS")
pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 30, MatrixOfVelocity=velocity)
tv1$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv1$origin$calculatedPseudotimeNoRoot<-which.min(tv1$pseudotime$calculatedPseudotimeNoRoot$res)
Walks(tv1,N=1000,origin_name = "calculatedPseudotimeNoRoot")

# pseudotime_50 <- get_pseudotime_from_velocity(tv1, 50, MatrixOfVelocity=velocity, tv1$origin$`A_B-blasts`)
# tv1$pseudotime$calculatedPseudotime50neighbours$res<-as.numeric(pseudotime_50)
# tv1$origin$calculatedPseudotime50neighbours<-tv1$origin$`A_B-blasts`
# Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours")
# 
# pseudotime_50_1955 <- get_pseudotime_from_velocity(tv1, 50, MatrixOfVelocity=velocity, as.integer(1955))
# tv1$pseudotime$calculatedPseudotime50neighbours1955$res<-as.numeric(pseudotime_50_1955)
# tv1$origin$calculatedPseudotime50neighbours1955<-as.integer(1955)
# Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours1955")
singler_labels <- readRDS("10X_analysis/singleR_annotation.RDS")
tv1$labels<-singler_labels
launch_shiny(tv1)

test<- as.matrix(readRDS( "/data/10X_analysis/all_X_PCA.RDS"))
tv2<-tviblindi(data=as.matrix(readRDS( "/data/10X_analysis/all_X_PCA.RDS")),labels=rep("A_B_blasts",4619))
#tv2<-tviblindi(data=as.matrix(data)[both$seurat_order,],labels=group_id[both$seurat_order])
DimRed(tv2)
DimRed(tv2,method="umap")

Set_origin(tv2,label = "A_B_blasts",origin_name = "A_B_blasts")


KNN(tv2)
Som(tv2,xdim = 15,ydim = 15) #kmeans clustering by default - 15*15 clusters
Filtration(tv2) #default setting is too conservative, less simplices could be created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv2,weighted = TRUE,origin_name = "A_B_blasts") 
Walks(tv2,N=1000,origin_name = "A_B_blasts") 

#tv1$vae$save("/data/patient2_0.5_subsample/patient2.json", "/data/patient2_0.5_subsample/patient2.h5")
##Adding un-rooted RNA velocity method
velocity <- readRDS("/data/10X_analysis/10x_PCA_velocity.RDS")
pseudotime_no_root <- get_pseudotime_from_velocity(tv2, 30, MatrixOfVelocity=velocity)
tv2$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv2$origin$calculatedPseudotimeNoRoot<-which.min(tv2$pseudotime$calculatedPseudotimeNoRoot$res)
Walks(tv2,N=1000,origin_name = "calculatedPseudotimeNoRoot")


launch_shiny(tv2)
