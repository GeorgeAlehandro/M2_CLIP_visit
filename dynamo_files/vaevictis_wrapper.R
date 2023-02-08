vv = reticulate::import("vaevictis")
#model <- vv$loadModel("/data/results/models/config.json","/data/results/models/weights.h5")
model <- vv$loadModel("/data/patient2_0.5_subsample/patient2.json", "/data/patient2_0.5_subsample/patient2.h5")

#layout <- model[[2]](tv1$data)
#plot(model[[2]](tv1$data))
#groups<-readRDS("/data/results/Hematopoiesis/hematopoiesis_cell_type.RDS")
#velocity<-readRDS("/data/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")
groups<-group_id
#velocity<-readRDS("/data/dynamo_files/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")
plot(model[[2]](velocity))
originalPCA_plus_velocity <- tv1$data + 1.2*velocity 
vaevictis_points <- as.data.frame(model[[2]](tv1$data))
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
