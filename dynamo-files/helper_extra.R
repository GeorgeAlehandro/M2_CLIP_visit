create_graph_from_transition_matrix <- function(tv1, n) {
  transition_matrix <- transition.matrix(tv1, n) 
  G <- graph_from_adjacency_matrix(transition_matrix,weighted = TRUE,mode = "undirected")
  listOfEdges <- split(as_edgelist(G), seq(nrow(as_edgelist(G))))
  names(listOfEdges)<-NULL
  return(list(G, listOfEdges))
}

calculate_complex_and_boundary <- function(flt) {
  cmplx<-build_boundary_CuR(flt)
  BB<-complex_to_boundaryF(cmplx = cmplx)
  B<-BB[[1]]
  B
}

calculate_Lap_and_Astar <- function(B, D0, D1) {
  Lap<-D0%*%B%*%D1%*%t(B)
  Astar<-D0%*%B%*%D1
  return(list(Lap, Astar))
}

calculate_alpha <- function(Lap, Astar, VD, IndexOfOriginCell) {
  right_hand <- -(Astar%*%VD)[-IndexOfOriginCell]
  left_hand <- Lap[-IndexOfOriginCell,-IndexOfOriginCell]
  alpha <- solve(left_hand, right_hand)
  return(alpha)
}


get_pseudotime_from_velocity <- function(tv1,nearest_neighbour_number, IndexOfRootCell) {
  # create graph and extract edges from transition matrix
  graph <- create_graph_from_transition_matrix(tv1, nearest_neighbour_number)
  G <- graph[[1]]
  listOfEdges <- graph[[2]]
  # calculate Laplacian and RHS
  laplacian_result <- build_Laplacian_and_RHS(transition_matrix, IndexOfRootCell)
  
  # calculate complex and boundary
  flt<-c(1:nrow(tv1$data),listOfEdges)
  B<-calculate_complex_and_boundary(flt)
  
  # load data
  V<-readRDS("/data/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")
  X<-tv1$data
  VD<-deRahmMap1f(B = B,X = X,V = V)
  
  D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
  D1<-Matrix::Diagonal(x=E(G)$weight)
  
  # calculate Lap and A*
  Lap_and_Astar <- calculate_Lap_and_Astar(B, D0, D1)
  Lap <- Lap_and_Astar[[1]]
  Astar <- Lap_and_Astar[[2]]
  # calculate alpha
  alpha_without_origin <- calculate_alpha(Lap, Astar, VD, IndexOfRootCell)
  alpha <-rep(0,nrow(tv1$data))
  alpha[-tv1$origin$HSC_hitting_time]<-as.numeric(alpha_without_origin)
  return(alpha)
}



pseudotime_50 <- get_pseudotime_from_velocity(tv1, 50, tv1$origin$HSC_hitting_time)
pseudotime_10 <- get_pseudotime_from_velocity(tv1, 10, tv1$origin$HSC_hitting_time)
pseudotime_30 <- get_pseudotime_from_velocity(tv1, 30, tv1$origin$HSC_hitting_time)
pseudotime_90 <- get_pseudotime_from_velocity(tv1, 90, tv1$origin$HSC_hitting_time)
tv1$pseudotime$calculatedPseudotime50neighbours$res<-as.numeric(pseudotime_50)
tv1$origin$calculatedPseudotime50neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime10neighbours$res<-as.numeric(pseudotime_10)
tv1$origin$calculatedPseudotime10neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime30neighbours$res<-as.numeric(pseudotime_30)
tv1$origin$calculatedPseudotime30neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime90neighbours$res<-as.numeric(pseudotime_90)
tv1$origin$calculatedPseudotime90neighbours<-tv1$origin$HSC_hitting_time

Walks(tv1,N=1000,origin_name = "calculatedPseudotime30neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime10neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime90neighbours")
