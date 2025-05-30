#' Generates null models for tripartite network
#'
#' A wrapper function to generate different null models for binary and quantitative tripartite networks
#'
#' @param trinet An 'igraph' object that represents a tripartite network.
#' @param null_N The number of null models to be generated.  Default to 100.
#' @param null_type  Character. Should be one of "sauve","sub_P","sub_Q" and "both_sub".See details.
#' @param sub_method The method to shuffle subnetworks. Must be provided when null_type ="sub_P","sub_Q" or "both_sub" .  a character specifying the null model algorithm listed on the help page of vegan::commsim. If null_type = 'sauve', it will be ignored.

#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes, b-nodes, c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#'
#' \strong{null_type}
#' \itemize{
#' \item “sub_P”, “sub_Q” and "both_sub" use the null model algorithms from \strong{vegan::commsim} to shuffle single subnetwork or both of them.
#' \item “sauve” rearranges connector species without changing subnetworks, following sauve et al.(2016).
#' }
#'
#' @return
#'
#' Return a list of null models of a tripartite network.
#'
#' @references
#' Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
#'
#' @importFrom igraph V
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @srrstats {G2.3} For univariate character input.
#' @srrstats {G2.3a} It provide four options of parameter: null_type.
#' @export
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#'
#' set.seed(123)
#' tri_null_list<-tri_null(Net,null_type="both_sub",sub_method="quasiswap")
#' set.seed(123)
#' tri_null_list<-tri_null(Net,null_type="sauve")
tri_null<-function(trinet, null_N=100, null_type=c("sauve","sub_P","sub_Q","both_sub"), sub_method){
   if(!null_type%in%c("sauve","sub_P","sub_Q","both_sub")|length(null_type)!=1){stop("Wrong input for null_type")}
   if(null_type%in%c("sub_P","sub_Q","both_sub")&missing(sub_method)){stop("sub_method should be provided for subnetwork null models")}
   mat<-as.matrix(trinet[])
   matP<-mat[V(trinet)$level==0,V(trinet)$level==1]
   matQ<-mat[V(trinet)$level==1,V(trinet)$level==2]
   if(null_type=="sub_P"){
      nm<-vegan::nullmodel(matP,sub_method)
      null_list<-simulate(nm, nsim=null_N)
      MAT<-mat
      null_network<-apply(null_list,3, function(x){
         MAT[V(trinet)$level==0,V(trinet)$level==1]<-x
         Nnetwork<-graph_from_adjacency_matrix(MAT,mode="max")#max to create an undirected graph
         V(Nnetwork)$name<-V(trinet)$name
         V(Nnetwork)$level<-V(trinet)$level
         dd<-igraph::layout_with_sugiyama(Nnetwork,layers=igraph::V(Nnetwork)$level)$layout
         dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
         dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
         dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
         Nnetwork$layout<-dd
         return(Nnetwork)
      })
   }
   else if(null_type=="sub_Q"){
      nm<-vegan::nullmodel(matQ,sub_method)
      null_list<-simulate(nm, nsim=null_N)
      MAT<-mat
      null_network<-lapply(null_list,3, function(x){
         MAT[V(network)$level==1,V(network)$level==2]<-x
         Nnetwork<-graph_from_adjacency_matrix(MAT,mode="max")
         V(Nnetwork)$name<-V(network)$name
         V(Nnetwork)$level<-V(network)$level
         dd<-igraph::layout_with_sugiyama(Nnetwork,layers=igraph::V(Nnetwork)$level)$layout
         dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
         dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
         dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
         Nnetwork$layout<-dd
         return(Nnetwork)
      })
   }
   else if(null_type=="both_sub"){
      nm_P<-vegan::nullmodel(matP,sub_method)
      null_P_list<-simulate(nm_P, nsim=null_N)
      nm_Q<-vegan::nullmodel(matQ,sub_method)
      null_Q_list<-simulate(nm_Q, nsim=null_N)
      MAT<-mat
      null_network<-lapply(1:null_N,function(x){
         MAT[V(trinet)$level==0,V(trinet)$level==1]<-null_P_list[,,x]
         MAT[V(trinet)$level==1,V(trinet)$level==2]<-null_Q_list[,,x]
         Nnetwork<-graph_from_adjacency_matrix(MAT,mode="max")
         V(Nnetwork)$name<-V(trinet)$name
         V(Nnetwork)$level<-V(trinet)$level
         dd<-igraph::layout_with_sugiyama(Nnetwork,layers=igraph::V(Nnetwork)$level)$layout
         dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
         dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
         dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
         Nnetwork$layout<-dd
         return(Nnetwork)
      })
   }
   else if(null_type=="sauve"){
      null_list1<-SauveR(null_N,matP,type="col")
      null_list2<-SauveR(null_N,matQ,type="row")
      MAT<-mat
      null_network<-lapply(1:null_N,function(x){
         MAT[V(trinet)$level==0,V(trinet)$level==1]<-null_list1[[x]]
         MAT[V(trinet)$level==1,V(trinet)$level==2]<-null_list2[[x]]
         Nnetwork<-graph_from_adjacency_matrix(MAT,mode="max")
         V(Nnetwork)$name<-V(trinet)$name
         V(Nnetwork)$level<-V(trinet)$level
         dd<-igraph::layout_with_sugiyama(Nnetwork,layers=igraph::V(Nnetwork)$level)$layout
         dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
         dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
         dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
         Nnetwork$layout<-dd
         return(Nnetwork)
      })
   }
   else
      stop("Error: null_type is not a valid input!")

   return(null_network)
}
