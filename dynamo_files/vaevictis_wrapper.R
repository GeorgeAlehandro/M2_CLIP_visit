vv = reticulate::import("vaevictis")
#model <- vv$loadModel("/data/results/models/config.json","/data/results/models/weights.h5")
model <- vv$loadModel("/data/10X_analysis/model.json", "/data/10X_analysis/model.h5")

#layout <- model[[2]](tv1$data)
#plot(model[[2]](tv1$data))
#groups<-readRDS("/data/results/Hematopoiesis/hematopoiesis_cell_type.RDS")
#velocity<-readRDS("/data/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")
groups<-tv0$labels[[1]]
#velocity<-readRDS("/data/dynamo_files/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")
plot(model[[2]](velocity))
originalPCA_plus_velocity <- tv3$data + velocity 
vaevictis_points <- as.data.frame(model[[2]](tv3$data))
vaevictis_velocity <- as.data.frame(model[[2]](originalPCA_plus_velocity))
vaevictis_velocity$celltypes <- groups

library(ggplot2)
ggplot()+
  geom_point(data=vaevictis_velocity,aes(V1,V2,color=celltypes))

# Create a data frame with the start and end points of the arrows
data_arrows <- data.frame(x1 = vaevictis_points$V1, y1 = vaevictis_points$V2, x2 = vaevictis_velocity$V1, y2 = vaevictis_velocity$V2)

# Combine the plots and draw the arrows
t <-geom_segment(data = data_arrows, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, "cm")), size = .85)


ggplot() + geom_point(data=vaevictis_points,aes(V1,V2,color="lightgrey")) + t + 
  theme(legend.position = "none")
vaevictis_points<-as.matrix(vaevictis_points)
#If need to reload
#rownames(vaevictis_points) <- NULL
#tv1$layout$lay<-vaevictis_points
ggplot()+
  geom_point(data=as.data.frame(data1),aes(PC_1,PC_2,color=vaevictis_velocity$celltypes))
