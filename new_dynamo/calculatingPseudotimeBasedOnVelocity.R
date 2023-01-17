devtools::install_github("stuchly/HodgePaths", ref= "main")
library("HodgePaths")

transition_matrix <-transition.matrix(tv1, 10) 
G<-graph_from_adjacency_matrix(transition_matrix,weighted = TRUE,mode = "undirected")
listOfEdges <- split(as_edgelist(G), seq(nrow(as_edgelist(G))))
names(listOfEdges)<-NULL

laplacian_result <- build_Laplacian_and_RHS(transition_matrix, tv1$origin$HSC_hitting_time)

flt<-c(1:nrow(tv1$data),listOfEdges)
cmplx<-build_boundary_CuR(flt)
BB<-complex_to_boundaryF(cmplx = cmplx)
B<-BB[[1]]

V<-readRDS("/data/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")

X<-tv1$data
VD<-deRahmMap1f(B = B,X = X,V = V)


D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
D1<-Matrix::Diagonal(x=E(G)$weight)

Lap<-D0%*%B%*%D1%*%t(B)
laplacian_result$Lap - Lap[-tv1$origin$HSC_hitting_time,-tv1$origin$HSC_hitting_time] # ~ 0
A<-t(B)
Astar<-D0%*%B%*%D1
right_hand <- -(Astar%*%VD)[-tv1$origin$HSC_hitting_time]
left_hand <- Lap[-tv1$origin$HSC_hitting_time,-tv1$origin$HSC_hitting_time]

alpha <- solve(left_hand, right_hand)
pseudo<-rep(0,nrow(tv1$data))
pseudo[-tv1$origin$HSC_hitting_time]<-as.numeric(alpha)
tv1$pseudotime$new_pseudotime$res<-as.numeric(pseudo1)
tv1$origin$new_pseudotime<-tv1$origin$HSC_hitting_time

Walks(tv1,N=1000,origin_name = "new_pseudotime")
