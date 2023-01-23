##Depends on tviblindi
transition.matrix<-function (x, K = 30, kernel = "SEMer", sym = "max",kepsilon = NULL) 
{
  
  if (K > dim(x$KNN$IND)[2]) {
    K <- min(K, dim(x$KNN)[2])
    warning("K > dim(KNN)[2]; K<-min(K,dim(x$KNN)[2])")
  }
  d <- KofRawN(x$KNN, K)
  d <- knn.raw2adj(d)
  dsym <- knn.spadj2sym(knn.adj2spadj(d))
  symB <- TRUE
  if (sym == "none") {
    sim <- knn.adj2spadjsim(d, kernel = kernel, epsilon = kepsilon)
    symB <- FALSE
  }
  else if (sym == "mean") 
    sim <- knn.spadj.symmetrize(knn.adj2spadjsim(d, kernel = kernel, 
                                                 epsilon = kepsilon))
  else if (sym == "prob") 
    sim <- knn.spadj.symmetrize.P(knn.adj2spadjsim(d, kernel = kernel, 
                                                   epsilon = kepsilon))
  else if (sym == "max") 
    sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel, 
                                          epsilon = kepsilon))
  else if (sym == "min") {
    d <- t(summary(dsym))
    sim <- knn.adj2spadjsim1(d, kernel = kernel, epsilon = kepsilon)
    symB = FALSE
  }
  else stop("symmetrisation not implemented")
  
  # .DD <- Matrix::rowSums(sim)
  # 
  # .DD <- Matrix::Diagonal(x=.DD^-1)  
  return(sim)
}

build_Laplacian_and_RHS<-function(A,origin,weights=NULL){
  .DD <- Matrix::rowSums(A)
  D <- Matrix::Diagonal(x = .DD)
  Dm <- Matrix::Diagonal(x = (.DD)^(-1))
  L <- Dm %*% (D - A)
  rm(D)
  unlabeled <- which(!(1:nrow(L) %in% origin))
  L <- L[unlabeled, unlabeled]
  if (is.null(weights)) {
    B <- matrix(1, nrow = length(unlabeled))
  }
  else {
    B <- Matrix::rowSums((Dm %*% 
                            A) * weights)
    B <- matrix(B[unlabeled], nrow = length(unlabeled))
  }
  return(list(Lap=L,RHS=B))
}

build_sym_Laplacian_and_RHS<-function(A,origin,weights=NULL){
  .DD <- Matrix::rowSums(A)
  D <- Matrix::Diagonal(x = .DD)
  Dm <- Matrix::Diagonal(x = (.DD)^(-1/2))
  L <- Dm %*% (D - A) %*% Dm
  rm(D)
  unlabeled <- which(!(1:nrow(L) %in% origin))
  L <- L[unlabeled, unlabeled]
  if (is.null(weights)) {
    B <- matrix(.DD[unlabeled]^(1/2), nrow = length(unlabeled))
  }
  else {
    B <- Matrix::rowSums((Matrix::Diagonal(x = (.DD)^(-1/2)) %*% 
                            A) * weights)
    B <- matrix(B[unlabeled], nrow = length(unlabeled))
  }
  return(list(Lap=L,RHS=B))
}
