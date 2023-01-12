
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
transition_matrix <- transition.matrix(tv1) 
D0<-Matrix::Diagonal(x=rowSums(transition_matrix))
D1<-Matrix::Diagonal(x=E(G)$weight)

#derahmap <- deRahmMap1f(get_boundary_matrix(tv1), tv1$data)


laplacian_result <- build_Laplacian_and_RHS(transition_matrix, tv1$origin$HSC_hitting_time)
# And weights?
laplacian_result <- build_Laplacian_and_RHS(BB, tv1$origin$HSC_hitting_time, weights = E(G)$weight)

laplacian_result$Lap * ALPHA = 

