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
  diag(sim)<-1
  
  .DD <- Matrix::rowSums(sim)
  
  .DD <- Matrix::Diagonal(x=.DD^-1)
  return(.DD%*%sim)
}