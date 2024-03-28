

#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph get.edgelist
#' @importFrom igraph V
#'
mEdges<-function(network){
   mat<-as.matrix(network[])
   #mat<-mat+t(mat)
   mat[mat>0]<-1
   net<-graph_from_adjacency_matrix(mat)
   EDGE1<-get.edgelist(net)
   EDGE1<-as.data.frame(EDGE1)
   EDGE1[,3]<-EDGE1[,2]
   EDGE1[,2]<-0
   EDGE1[,4]<-0
   EDGE1[,5]<-1
   EDGE1[EDGE1[,1]%in%V(network)$name[V(network)$level==0],2]<-1
   EDGE1[EDGE1[,1]%in%V(network)$name[V(network)$level==1],2]<-2
   EDGE1[EDGE1[,1]%in%V(network)$name[V(network)$level==2],2]<-3
   EDGE1[EDGE1[,3]%in%V(network)$name[V(network)$level==0],4]<-1
   EDGE1[EDGE1[,3]%in%V(network)$name[V(network)$level==1],4]<-2
   EDGE1[EDGE1[,3]%in%V(network)$name[V(network)$level==2],4]<-3
   for (i in 1:length(V(network)$name)) {
      EDGE1[EDGE1[,1]%in%V(network)$name[i],1]<-i
      EDGE1[EDGE1[,3]%in%V(network)$name[i],3]<-i
   }
   EDGE1[,1]<-as.numeric(EDGE1[,1])
   EDGE1[,3]<-as.numeric(EDGE1[,3])
   return((EDGE1))
}
