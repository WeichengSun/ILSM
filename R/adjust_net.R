#' Adjust the tripartite network by deleting non-connector nodes in shared set
#'
#' @param network The original tripartite network
#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#' The original tripartite network should have a 'level' attribute for nodes with 0 for a-nodes, 1 for b-nodes and 2 for c-nodes.
#'
#' @importFrom igraph V
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @noRd
#'
#' @srrstats {G1.4a} This is a non-exported function.
#'
#' @return
#'
#' An 'igraph' object representing tripartite network that contains only connector nodes/species in the shared set.

adjust_net<-function(network,weighted=T){
  PHP1<-as.matrix(network[])
  PH1<-PHP1[(V(network)$level)==0,(V(network)$level)==1]
  HP1<-PHP1[(V(network)$level)==1,(V(network)$level)==2]
  M111_1<-(apply(abs(PH1),2,sum)*apply(abs(HP1),1,sum))
  if(sum(M111_1>0)){
    PHP_test<-PHP1[!rownames(PHP1)%in%names(M111_1)[(M111_1==0)],!rownames(PHP1)%in%names(M111_1)[(M111_1==0)]]
    # apply(PHP_test[(rowSums(PH)),],1,function(x){})
    PHP_test1<-PHP_test[!rownames(PHP_test)%in%rownames(PH1)[rowSums(PHP_test[rownames(PH1),])==0],!rownames(PHP_test)%in%rownames(PH1)[rowSums(PHP_test[rownames(PH1),])==0]]
    PHP_test2<-PHP_test1[!rownames(PHP_test1)%in%colnames(HP1)[colSums(PHP_test1[,colnames(HP1)])==0],!rownames(PHP_test1)%in%colnames(HP1)[colSums(PHP_test1[,colnames(HP1)])==0]]
    NODE_logi<-(!rownames(PHP1)%in%c(names(M111_1)[(M111_1==0)],rownames(PH1)[rowSums(PHP_test[rownames(PH1),])==0],colnames(HP1)[colSums(PHP_test1[,colnames(HP1)])==0]))
    PHP_test2<-igraph::graph_from_adjacency_matrix(PHP_test2,weighted = weighted)
    V(PHP_test2)$name<-(V(network)$name)[NODE_logi]
    V(PHP_test2)$level<-(V(network)$level)[NODE_logi]
    return(PHP_test2)
  }
  else
    stop("No connector nodes.")
}
