#' Generating an example tripartite network randomly
#'
#' Generate a binary tripartite network with three groups of nodes (a-, b- and c-nodes) and two subnetworks (P and Q). Subnetwork P contains links between a- and b-nodes; Subnetwork Q contains links between b- and c-nodes;b-nodes (some of which are connector nodes) are shared nodes between two subnetworks.
#'
#' @param N_a  The number of nodes in the a-node group.
#' @param N_b  The number of nodes in the b-node group.
#' @param N_c  The number of nodes in the c-node group.
#' @param Co  The probability of creating a link between any two nodes. It ranges from 0 to 1.
#' @param output_matrices Logical. Whether to output the entire adjacency matrix of the network and subnetworks. Defaults to FALSE.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph V
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout_with_sugiyama
#' @importFrom stats runif
#'
#' @return
#' Return a random binary tripartite network.
#' @export
#'
#' @srrstats {G1.0} This standard belongs here.
#' @srrstats {G2.0} It can input single- or multi- data.
#' @srrstats {G2.0a} It provide explicit secondary documentation.
#' @srrstats {G2.4,G2.4b} This standard belongs here.
#' @srrstats {G2.6} It accepts one-dimensional input.
#'
#' @references
#'
#' Pilosof, S., Porter, M., Pascual, M. et al. The multilayer nature of ecological networks. Nat Ecol Evol 1, 0101 (2017). https://doi.org/10.1038/s41559-017-0101
#'
#' @examples
#'
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#' plot(Net)
#'
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2,output_matrices=T)
#' Net
#'
build_toy_net<-function(N_a, N_b, N_c, Co,output_matrices=FALSE){
   if(N_a<3||N_b<3||N_c<3)
      stop("Error: please make sure N_a>=3, N_b>=3 and N_c>=3!!!")
   lay<-N_a+N_b+N_c
   network<-matrix(0,lay,lay)
   node<-1:lay
   network[1:N_a,(N_a+1):(N_a+N_b)]<-as.numeric(runif(N_a*N_b)<Co)
   network[(N_a+1):(N_a+N_b),(N_a+N_b+1):(N_a+N_b+N_c)]<-as.numeric(runif(N_b*N_c)<Co)
   network<-igraph::graph_from_adjacency_matrix(network,"max")
   igraph::V(network)$name<-c(paste0("a",1:N_a),paste0("b",1:N_b),paste0("c",1:N_c))
   igraph::V(network)$level<-c(rep(0,N_a),rep(1,N_b),rep(2,N_c))
   dd<-igraph::layout_with_sugiyama(network,layers=igraph::V(network)$level)$layout
   dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
   dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
   dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
   network$layout<-dd
   if(output_matrices){
      if(inherits(network,"igraph")==F){
         warning("Warning: the network is not an igraph object!!!")
         return(network)
      }
      tri_mat<-as.matrix(network[])
      subnet_P<-tri_mat[(igraph::V(network)$level)==0,(igraph::V(network)$level)==1]
      subnet_Q<-tri_mat[(igraph::V(network)$level)==1,(igraph::V(network)$level)==2]
      return(list(network=network,supraadjacency_matrix= tri_mat,subnet_P=t(subnet_P),subnet_Q=subnet_Q))
   }
   return(network)
}
