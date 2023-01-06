library("Matrix")

shape <- c(1,2,3,4,5,6)
possible.columns = c()
selecting <- c(1,2,3)
#vertice <- setRefClass("vertice", fields=list(name="character", edges="vector"))
vertice <- setRefClass("vertice", fields=list(name="character", edges="list"))

addEdge<-function(object, node) {
  object$edges <- append(object$edges,node)
  return(object)
}

check<-function(object, node) {
  object$edges <- append(object$edges,node)
  return(object) 
}

node1 <- vertice$new(name="1")
node2 <- vertice$new(name="2")
node3 <- vertice$new(name="3")
node4 <- vertice$new(name="4")
node5 <- vertice$new(name="5")
node6 <- vertice$new(name="6")
node7 <- vertice$new(name="7")

node1<-addEdge(node1,node2)
node2<-addEdge(node2,node1)

node1<-addEdge(node1,node3)
node3<-addEdge(node3,node1)

node2<-addEdge(node2,node6)
node6<-addEdge(node6,node2)

node2<-addEdge(node2,node3)
node3<-addEdge(node3,node2)

node2<-addEdge(node2,node4)
node4<-addEdge(node4,node2)

node3<-addEdge(node3,node4)
node4<-addEdge(node4,node3)

node4<-addEdge(node4,node5)
node5<-addEdge(node5,node4)

node4<-addEdge(node4,node7)
node7<-addEdge(node7,node4)


node5<-addEdge(node5,node7)
node7<-addEdge(node7,node5)

node6<-addEdge(node6,node5)
node5<-addEdge(node5,node6)


containsNode<- function(node, listOfNodes){
  for (nodeInList in listOfNodes){
    if (node$name == nodeInList$name){
      return(TRUE)
    }
  }
    return(FALSE)
}
# for (i in 1:(length(shape)-1)){
#   
#   for (j in 1:length(selecting)){
#     
#   }
#   print(paste(shape[i],shape[i+1]))
#   
#   
#   possible.columns<- append(possible.columns, paste(shape[i],shape[i+1]))
# }
listOfNodes <- list(node1,node2, node3, node4, node5, node6, node7)
extractPairs <- function(node, returnedPairs){
  for (neighbourNode in node$edges){
    if (!(paste0("[",neighbourNode$name,", ",node$name,"]") %in% returnedPairs)){
    returnedPairs <- append(returnedPairs,paste0("[",node$name,", ",neighbourNode$name,"]"))
    }
  }
  return(returnedPairs)
}
extractedPairs <- c()
for (eachNode in listOfNodes){
  extractedPairs <- extractPairs(eachNode, extractedPairs)
}


extractPairsFromTriangle <- function(node, returnedPairs, listOfNodesOfTriangle){
  for (neighbourNode in node$edges){
    if (containsNode(neighbourNode,listOfNodesOfTriangle)){
    if (!(paste0("[",neighbourNode$name,", ",node$name,"]") %in% returnedPairs)){
      returnedPairs <- append(returnedPairs,paste0("[",node$name,", ",neighbourNode$name,"]"))
    }
    }
  }
  return(returnedPairs)
}


triangle123 <- list(node1, node2, node3)
extractedPairsFromTriangle <- c()
for (eachNodeInTriangle in triangle123){
  extractedPairsFromTriangle <- extractPairsFromTriangle(eachNodeInTriangle, extractedPairsFromTriangle,triangle123)
}

matchVectors(pairs, query){
  
}
nameddd <- as.vector("[1,2,3]")
df <- data.frame(vectorww = extractedPairs, Shape123 = extractedPairs %in% extractedPairsFromTriangle)
