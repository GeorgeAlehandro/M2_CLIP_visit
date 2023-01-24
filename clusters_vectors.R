max(tv1$clusters)#
plot(tv1$layout$lay[,1], tv1$layout$lay[,2])
mean_clusters <- data.frame()
library(dplyr)
# After running vaevictis_wrapper.R
data <- data.frame(dim_red = tv1$layout$lay, clusters = tv1$clusters, velocity = vaevictis_velocity)

mean_cluster <- data %>%
  group_by(clusters) %>%
  summarise_all(funs(mean))

data_arrows <- data.frame(x1 = mean_cluster[,2], y1 = mean_cluster[,3], x2 = mean_cluster[,4], y2 = mean_cluster[,5])

t <-geom_segment(data = data_arrows, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, "cm")), size = .85)


ggplot() + geom_point(data=mean_cluster,aes(dim_red.V1,dim_red.V2,color="lightgrey"), size = 5.5) + t + 
  theme(legend.position = "none")


mean_cluster <- as.data.frame(mean_cluster)
plot(mean_cluster[,2], mean_cluster[,3])
