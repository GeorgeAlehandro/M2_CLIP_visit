get_boundary_matrix <- function(tviblindi_object){
  A<-transition.matrix(tviblindi_object)
  G<-graph_from_adjacency_matrix(A,weighted = TRUE,mode = "undirected")
  xy.list <- split(as_edgelist(G), seq(nrow(as_edgelist(G))))
  names(xy.list)<-NULL
  flt<-c(1:nrow(tv1$data),xy.list)
  cmplx<-build_boundary_CuR(flt)
  BB<-complex_to_boundaryF(cmplx = cmplx)
  return (BB)
}


get_weights <- function(tviblindi_object){
  A<-transition.matrix(tviblindi_object)
  G<-graph_from_adjacency_matrix(A,weighted = TRUE,mode = "undirected")
  return(E(G)$weight)
}

get_pseudotime_from_velocity <- function(tviblindi_object, nearest_neighbour_number){
  
  transition_matrix <-transition.matrix(tviblindi_object, nearest_neighbour_number) 
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
  right_hand <- (Astar%*%VD)[-tv1$origin$HSC_hitting_time]
  left_hand <- Lap[-tv1$origin$HSC_hitting_time,-tv1$origin$HSC_hitting_time]
  
  alpha <- solve(left_hand, right_hand)

  return(alpha)
}
pseudotime_50 <- get_pseudotime_from_velocity(tv1, 50)
pseudotime_10 <- get_pseudotime_from_velocity(tv1, 10)
pseudotime_30 <- get_pseudotime_from_velocity(tv1, 30)
pseudotime_90 <- get_pseudotime_from_velocity(tv1, 90)
tv1$pseudotime$calculatedPseudotime50neighbours$res<--as.numeric(pseudotime_50)
tv1$origin$calculatedPseudotime50neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime10neighbours$res<--as.numeric(pseudotime_10)
tv1$origin$calculatedPseudotime10neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime30neighbours$res<--as.numeric(pseudotime_30)
tv1$origin$calculatedPseudotime30neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime90neighbours$res<--as.numeric(pseudotime_90)
tv1$origin$calculatedPseudotime90neighbours<-tv1$origin$HSC_hitting_time

Walks(tv1,N=1000,origin_name = "calculatedPseudotime30neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime10neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime90neighbours")
