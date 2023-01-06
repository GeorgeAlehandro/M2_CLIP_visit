library("Matrix")

shape <- c(1,2,3,4,5,6)
possible.columns = c()
selecting <- c(1,2,3)

setClass("node", slots=list(name="character",
                            edges="list"))

# Calling object
node

addEdge<-           function(object, node) {
  object@edges <- append(object@edges,node)
  return(object) 
}
setMethod("show",
          "node",
          function(object) {
            print(object@name)
          }
)
setMethod("addEdge",
          "node",
          function(object, node) {
            object@edges <- append(object@edges,node)
            print("test")
            print(object@edges)
            return(object)           
            }
)
obj1 <- new("node",name="1")
obj2 <- new("node",name="2")
obj3 <- new("node",name="3")
obj4 <- new("node",name="4")
obj5 <- new("node",name="5")
obj6 <- new("node",name="6")
obj7 <- new("node",name="7")

obj1<-addEdge(obj1,obj2)
obj2<-addEdge(obj2,obj1)
  
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
