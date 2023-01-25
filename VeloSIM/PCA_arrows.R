dataset<-readRDS("/data/VeloSIM/datasets/2000_100_0.2.rds")
data_counts <- dataset$counts_s_PCA[,1:2]
velocity<-dataset$velocity_100_pca[,1:2]
groups<-dataset$backbone
plot(velocity[,1], velocity[,2])
originalPCA_plus_velocity <- data_counts + velocity 
vaevictis_points <- as.data.frame(data_counts)
vaevictis_velocity <- as.data.frame(originalPCA_plus_velocity)
vaevictis_velocity$celltypes <- groups
vaevictis_points$celltypes <- groups
library(ggplot2)
ggplot()+
  geom_point(data=vaevictis_points,aes(PC1,PC2,color=celltypes))

# Create a data frame with the start and end points of the arrows
data_arrows <- data.frame(x1 = vaevictis_points$PC1, y1 = vaevictis_points$PC2, x2 = vaevictis_velocity$PC1, y2 = vaevictis_velocity$PC2)

# Combine the plots and draw the arrows
t <-geom_segment(data = data_arrows, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, "cm")), size = .85)


ggplot() + geom_point(data=vaevictis_points,aes(PC1,PC2,color="lightgrey")) + t + 
  theme(legend.position = "none")
vaevictis_points<-as.matrix(vaevictis_points)
rownames(vaevictis_points) <- NULL
tv1$layout$lay<-vaevictis_points
