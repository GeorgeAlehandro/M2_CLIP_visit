
devtools::install_github("stuchly/HodgePaths",ref="main")
library("HodgePaths")
create_graph_from_transition_matrix <- function(tv1, n) {
  transition_matrix <- transition.matrix(tv1, n) 
  G <- graph_from_adjacency_matrix(transition_matrix,weighted = TRUE,mode = "undirected")
  listOfEdges <- split(as_edgelist(G), seq(nrow(as_edgelist(G))))
  names(listOfEdges)<-NULL
  return(list(transition_matrix,G, listOfEdges))
}

calculate_complex_and_boundary <- function(flt) {
  cmplx<-build_boundary_CuR(flt)
  BB<-complex_to_boundaryF(cmplx = cmplx)
  B<-BB[[1]]
  B
}
# Intermediate 
calculate_standard_Astar <- function(B) {
  Astar<-B
  return(Astar)
}
calculate_standard_Lap <- function(B) {
  Lap<-B%*%t(B)
  return(Lap)
}

calculate_standard_Lap_sym <- function(B, D0, D1) {
  Dm <- Matrix::Diagonal(x=diag(D0)^(1/2))
  Lap_sym<-Dm%*%B%*%D1%*%t(B)%*%Dm
  # Lap_sym <- Dm%*%Lap%*%Dm_minus1
  return(Lap_sym)
}

calculate_alpha <- function(Lap, Astar, VD, IndexOfOriginCell,lap_symmetrical) {
  right_hand <- -(Astar%*%VD)[-IndexOfOriginCell]
  left_hand <- Lap[-IndexOfOriginCell,-IndexOfOriginCell]
  alpha <- solve(left_hand, right_hand)
  return(alpha)
}

calculate_alpha_pseudo_inverse <- function(Lap, Astar, D0, VD,lap_symmetrical=F) {
  if(!lap_symmetrical){
    right_hand <- -(Astar%*%VD)
    left_hand <- Lap
    alpha<-HodgePaths:::bicgSparse(left_hand,as.numeric(right_hand),nb_iter = 1500)
    alpha<-alpha$x
  }
  else{
    Dm_minus <- Matrix::Diagonal(x=diag(D0)^(-1/2))
    Dm_minus_inverse <- Matrix::Diagonal(x=diag(Dm_minus)^(-1))
    left_hand <<- Lap
    right_hand <<- -(Dm_minus%*%Astar%*%VD)
    epsilon<<-HodgePaths:::bicgSparse(left_hand,as.numeric(right_hand), nb_iter = 1500)
    alpha <- Dm_minus_inverse%*%epsilon$x
    #  epsilon<<-solve(left_hand,as.numeric(right_hand))
    # alpha <- Dm_minus_inverse%*%epsilon
  }
  return(alpha)
}



get_pseudotime_from_velocity_standard_lap <- function(tv1, nearest_neighbour_number=30, IndexOfRootCell=NULL, MatrixOfVelocity = "RDS", lap_symmetrical = F) {
  # create graph and extract edges from transition matrix
  graph <- create_graph_from_transition_matrix(tv1, nearest_neighbour_number)
  transition_matrix <- graph[[1]]
  G <<- graph[[2]]
  listOfEdges <- graph[[3]]
  # calculate complex and boundary
  flt<-c(1:nrow(tv1$data),listOfEdges)
  B<<-calculate_complex_and_boundary(flt)
  
  # load data
  if (MatrixOfVelocity == "RDS"){
    V<-readRDS("/data/dynamo_files/results/Hematopoiesis/hematopoiesis_PCA_velocity.RDS")}
  else{
    V <- MatrixOfVelocity
  }
  X<-tv1$data
  VD<-deRahmMap1f(B = B,X = X,V = V)
  Astar <- B
  # calculate Lap and A*
  if (lap_symmetrical){
    Lap <- calculate_standard_Lap_sym(B)
  }
  else{
    Lap <<- calculate_standard_Lap(B)
  }
  # calculate alpha depending if the root cell is selected or not
  if (is.null(IndexOfRootCell)){
    print(lap_symmetrical)
    alpha <- calculate_alpha_pseudo_inverse(Lap, Astar, VD, D0=D0, lap_symmetrical=lap_symmetrical)
    return(alpha)
    
  }
  else{
    alpha_without_origin <- calculate_alpha(Lap, Astar, VD, IndexOfRootCell,lap_symmetrical=lap_symmetrical)
    alpha <-rep(0,nrow(tv1$data))
    alpha[-IndexOfRootCell]<-as.numeric(alpha_without_origin)
    return(alpha)
  }
  
}
pseudotime_30_standard <- get_pseudotime_from_velocity_standard_lap(tv1, 30, IndexOfRootCell=tv1$origin$`5_6_hitting_time`,MatrixOfVelocity = simulated_data$velocity_100_pca)
tv1$pseudotime$calculatedPseudotimeStandard$res<-as.numeric(pseudotime_30_standard)
tv1$origin$calculatedPseudotimeStandard<-tv1$origin$`5_6_hitting_time`
Walks(tv1,N=1000,origin_name = "calculatedPseudotimeStandard")
