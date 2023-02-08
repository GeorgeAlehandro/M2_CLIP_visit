
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
calculate_Astar <- function(B, D0, D1) {
  Astar<-D0%*%B%*%D1
  return(Astar)
}
calculate_Lap <- function(B, D0, D1) {
  Lap<-D0%*%B%*%D1%*%t(B)
  return(Lap)
}

calculate_Lap_sym <- function(B, D0, D1) {
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



get_pseudotime_from_velocity <- function(tv1, nearest_neighbour_number=30, IndexOfRootCell=NULL, MatrixOfVelocity = "RDS", lap_symmetrical = F) {
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
  D0<<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
  D1<<-Matrix::Diagonal(x=E(G)$weight)
  Astar <- calculate_Astar(B,D0,D1)
  # calculate Lap and A*
  if (lap_symmetrical){
    Lap <- calculate_Lap_sym(B, D0, D1)
  }
  else{
  Lap <- calculate_Lap(B, D0, D1)
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

pseudotime_no_root <- get_pseudotime_from_velocity(tv1, 30)
pseudotime_50 <- get_pseudotime_from_velocity(tv1, 50, tv1$origin$HSC_hitting_time)
pseudotime_50_sym <- get_pseudotime_from_velocity(tv1, 50, tv1$origin$HSC_hitting_time, lap_symmetrical = T)
pseudotime_10 <- get_pseudotime_from_velocity(tv1, 10, tv1$origin$HSC_hitting_time)
pseudotime_30 <- get_pseudotime_from_velocity(tv1, 30, tv1$origin$HSC_hitting_time)
pseudotime_90 <- get_pseudotime_from_velocity(tv1, 90, tv1$origin$HSC_hitting_time)

test <- get_pseudotime_from_velocity(tv1, 90, lap_symmetrical = T )
test_non_sym <- get_pseudotime_from_velocity(tv1, 90, lap_symmetrical = F )

tv1$pseudotime$symNoRoot$res<-as.numeric(test)
tv1$origin$symNoRoot<-which.min(tv1$pseudotime$symNoRoot$res)

tv1$pseudotime$NoRoot$res<-as.numeric(test_non_sym)
tv1$origin$NoRoot<-which.min(tv1$pseudotime$NoRoot$res)

tv1$pseudotime$calculatedPseudotimeNoRoot$res<-as.numeric(pseudotime_no_root)
tv1$origin$calculatedPseudotimeNoRoot<-which.min(tv1$pseudotime$calculatedPseudotimeNoRoot$res)



tv1$pseudotime$calculatedPseudotime50neighbours$res<-as.numeric(pseudotime_50)
tv1$origin$calculatedPseudotime50neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime50_SYM_neighbours$res<-as.numeric(pseudotime_50_sym)
tv1$origin$calculatedPseudotime50_SYM_neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime10neighbours$res<-as.numeric(pseudotime_10)
tv1$origin$calculatedPseudotime10neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime30neighbours$res<-as.numeric(pseudotime_30)
tv1$origin$calculatedPseudotime30neighbours<-tv1$origin$HSC_hitting_time

tv1$pseudotime$calculatedPseudotime90neighbours$res<-as.numeric(pseudotime_90)
tv1$origin$calculatedPseudotime90neighbours<-tv1$origin$HSC_hitting_time


Walks(tv1,N=1000,origin_name = "calculatedPseudotimeNoRoot")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime30neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime10neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime50neighbours")
Walks(tv1,N=1000,origin_name = "calculatedPseudotime90neighbours")


from_bound <- calculate_Lap(B,D0,D1)

from_lap <- build_Laplacian_and_RHS(transition.matrix(tv1), origin=tv1$origin$`5_6_hitting_time`)

from_bound[-tv1$origin$HSC_hitting_time, -tv1$origin$HSC_hitting_time] - from_lap[[1]]

from_bound <- calculate_Lap_sym(B,D0,D1)
from_lap <- build_sym_Laplacian_and_RHS(transition.matrix(tv1), origin=tv1$origin$HSC_hitting_time)
from_bound[-tv1$origin$HSC_hitting_time, -tv1$origin$HSC_hitting_time] - from_lap[[1]]
