
# library(tviblindi)
# data<-readRDS("/data/results/sc-velo-hematopoiesis/scVelo_hemato_pca.RDS")
# group_id<-readRDS("/data/results/Hematopoiesis/hematopoiesis_cell_type.RDS")
# tvvelo<-tviblindi(data=data,labels=group_id)
# DimRed(tvvelo)

vv = reticulate::import("vaevictis")
model <- vv$loadModel("/data/results/models/config_scvelo.json","/data/results/models/weights_scvelo.h5")
layout <- model[[2]](tv1$data)
plot(model[[2]](tv1$data))
groups<-readRDS("/data/results/Hematopoiesis/hematopoiesis_cell_type.RDS")
velocity<-readRDS("/data/results/sc-velo-hematopoiesis/scVelo_hemato_pca.RDS")*10^-1
plot(model[[2]](velocity))
originalPCA_plus_velocity <- tvvelo$data + velocity
vaevictis_points <- as.data.frame(model[[2]](tvvelo$data))
vaevictis_velocity <- as.data.frame(model[[2]](originalPCA_plus_velocity))
vaevictis_points$celltypes <- groups
vaevictis_velocity$celltypes <- groups

library(ggplot2)
ggplot()+
  geom_point(data=vaevictis_velocity,aes(V1,V2,color=celltypes))

ggplot()+
  geom_point(data=vaevictis_points,aes(V1,V2,color=celltypes))
p <- ggplot() + geom_point(data=vaevictis_points,aes(V1,V2,color=celltypes)) + t

# Create a data frame with the start and end points of the arrows
data_arrows <- data.frame(x1 = vaevictis_points$V1, y1 = vaevictis_points$V2, x2 = vaevictis_velocity$V1, y2 = vaevictis_velocity$V2)
# Combine the plots and draw the arrows
t <-geom_segment(data = data_arrows, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, "cm")))


p <- ggplot() + geom_point(data=vaevictis_points,aes(V1,V2,color=celltypes)) + t
