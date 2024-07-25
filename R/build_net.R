#' Generating tripartite network
#'
#' Generating a network of three layers. All layers of network contain \code{lay_0}, \code{lay_1} and \code{lay_2} nodes respectively.
#'
#' @param lay_0  The number of nodes in the first layer.
#' @param lay_1  The number of nodes in the second layer.
#' @param lay_2  The number of nodes in the third layer.
#' @param C_lay  The probability of each node interact with the other one. It ranges from 0 to 1.
#' @param asmatrices Logical. whether to output the overall adjacency matrix of the network and the corresponding interaction matrix of the respective subnetworks. Defaults to FALSE.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph V
#' @importFrom igraph plot.igraph
#' @importFrom igraph layout_with_sugiyama
#' @importFrom stats runif
#'
#' @return
#' Return a tripatite network of direction. The network contains three groups of species and interactions within layers, and there is no link among each group of nodes within one layer.
#' @export
#' @examples
#'
#' set.seed(12)
#' d <- build_net(11,15,16,0.2)
#' plot(d)
#'
#' set.seed(12)
#' N <- build_net(11,15,16,0.2,asmatrices=FALSE)
#' N
#'
build_net<-function(lay_0, lay_1, lay_2, C_lay, asmatrices=FALSE){
   if(lay_0<3||lay_1<3||lay_2<3)
      stop("Error: please make lay_0>=3, lay_1>=3 and lay_2>=3!!!")
   lay<-lay_0+lay_1+lay_2
   network<-matrix(0,lay,lay)
   node<-sample(lay)
   network[1:lay_0,(lay_0+1):(lay_0+lay_1)]<-as.numeric(runif(lay_0*lay_1)<C_lay)
   network[(lay_0+1):(lay_0+lay_1),(lay_0+lay_1+1):(lay_0+lay_1+lay_2)]<-as.numeric(runif(lay_1*lay_2)<C_lay)
   network<-igraph::graph_from_adjacency_matrix(network)
   igraph::V(network)$name<-paste0("S",node,collapse =NULL,recycle0 = F)
   igraph::V(network)$level<-c(rep(0,lay_0),rep(1,lay_1),rep(2,lay_2))
   dd<-igraph::layout_with_sugiyama(network,layers=V(network)$level)$layout
   dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
   dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
   dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
   network$layout<-dd
   if(asmatrices){
      if(inherits(network,"igraph")==F){
         warning("Warning: the type of network is not an igraph!!!")
         return(network)
      }
      PHP<-as.matrix(network[])
      PH<-PHP[(V(network)$level)==0,(V(network)$level)==1]
      HP<-PHP[(V(network)$level)==1,(V(network)$level)==2]
      return(list(network=network,supraadjacency_matrix=PHP,subnetwork1=t(PH),subnetwork2=HP))
   }
   return(network)
}
