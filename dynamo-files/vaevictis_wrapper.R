vv = reticulate::import("vaevictis")
model <- vv$loadModel("config.json","weights.h5")
layout <- model[[2]](tv1$data)
plot(model[[2]](tv1$data))
groups<-readRDS("/data/hematopoiesis_cell_type.RDS")
velocity<-readRDS("/data/hematopoiesis_PCA_velocity.RDS")
plot(model[[2]](velocity))
originalPCA_plus_velocity <- tv1$data + velocity
vaevictis_points <- as.data.frame(model[[2]](tv1$data))
vaevictis_velocity <- as.data.frame(model[[2]](originalPCA_plus_velocity))
vaevictis_points$celltypes <- groups
vaevictis_velocity$celltypes <- groups

line_data <- data.frame()
for(i in 1:nrow(vaevictis_points)) {
  # Insert a row from df1
  line_data <- rbind(line_data, vaevictis_points[i,])
  
  # Insert a row from df2
  line_data <- rbind(line_data, vaevictis_velocity[i,])
}
line_data

library(ggplot2)
ggplot()+
  geom_point(data=vaevictis_velocity,aes(V1,V2,color=celltypes))

ggplot()+
  geom_point(data=vaevictis_points,aes(V1,V2,color=celltypes))
p <- ggplot() 
for(i in seq(from=1, to=2000-1, by=2)){
  print(i)
  p <- p + geom_path(data = line_data[i:(i+1),], aes(V1, V2), linetype = "dotted")
}
p
ggplot()+
  geom_point(data=vaevictis_points,aes(V1,V2,color=celltypes))+
  p



