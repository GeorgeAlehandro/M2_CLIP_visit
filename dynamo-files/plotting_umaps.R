library(tidyverse)
library(umap)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library("ggplot2")
umap_coordinates <- as.data.frame(readRDS("/home/georgealehandro/dynamo-files/hematopoiesis_UMAP.RDS"))
colnames(umap_coordinates) <- c("UMAP1", "UMAP2")
groups <- readRDS("/home/georgealehandro/dynamo-files/hematopoiesis_cell_type.RDS")
umap_coordinates$celltypes <- groups
umap_coordinates %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = celltypes))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")


df_velocity <- as.data.frame(readRDS("/home/georgealehandro/dynamo-files/hematopoiesis_PCA_velocity.RDS"))


umap_velocity <- df_velocity %>%
  umap()

umap_velocity <- umap_velocity$layout
umap_velocity<- as.data.frame(umap_velocity)
colnames(umap_velocity) <- c("UMAP1", "UMAP2")
umap_velocity$celltypes <- groups
umap_velocity %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = celltypes))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP Velocity plot")

library(ggplot2)
line_data <- data.frame(x = c(umap_coordinates$UMAP1,umap_velocity$UMAP1), y = c(umap_coordinates$UMAP2,umap_velocity$UMAP2))
line_data <- data.frame()
for(i in 1:nrow(umap_coordinates)) {
  # Insert a row from df1
  line_data <- rbind(line_data, umap_coordinates[i,])
  
  # Insert a row from df2
  line_data <- rbind(line_data, umap_velocity[i,])
}
line_data



ggplot()+
  geom_point(data=umap_coordinates,aes(UMAP1,UMAP2,color=celltypes))+
  geom_point(data=umap_velocity,aes(UMAP1,UMAP2,color=celltypes))+ 
  geom_path(data = line_data, aes(UMAP1, UMAP2), linetype = "dotted", linewidth =0.5)

library(ggplot2)

# Create the data for the two plots
df1 <- data.frame(x = c(1, 2, 3), y = c(3, 4,5))

# Create a plot with a for loop
p <- ggplot()
for(i in 1:(nrow(df1)-1)) {
  p <- p + geom_path(data = df1[i:(i+1),], aes(x, y), linetype = "dotted")
}
p

