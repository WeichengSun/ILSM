#' Null model of multilayer network
#'
#' The null model could be generated according to different matrix scrambling algorithms for interconnection patterns in the multilayer network
#'
#' @param network A multilayer(tripartite) network of 'igraph' class. The network contains interlayer links and without intralayer links.
#' @param number A numeric value. The number of null model.  Default to NULL representing number 1.
#' @param null_type Logical. Four matrix scrambling algorithms. If null_type = NULL, default to "all".
#'
#' @details
#'
#' \strong{null_type}
#'
#' \itemize{
#' \item{For each of the four types of null models, there are corresponding algorithms. The first type, “subnetwork1”, involved scrambling the adjacency matrix of the first and second layers of the multilayer network.}
#' \item{The second type, “subnetwork2”, focused on scrambling the adjacency matrix of the second and third layers. }
#' \item{Comprehensively, the third type, “all”, blended the approaches of the first two to disarrange the entire network's adjacency matrix, achieving a thorough perturbation of the network's structure. }
#' \item{The last type named “Savue” that disarranged inherent structure in terms of the groups of species connected by each interconnecting species of every subnetworks, thus exhibiting different interconnection patterns.}}
#'
#' \strong{network}
#'
#' About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}
#'
#' @return
#'
#' Return a list contains one or more elements. Each element represent a null model of multilayer network.
#'
#' @references
#' Vázquez, D. P., C. J. Melian, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction strength in ecological networks. Oikos 116: 1120-1127.
#'
#' Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
#'
#' @importFrom igraph V
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @export
#'
#' @examples
#'
#' set.seed(12)
#' d <- build_net(11,22,21,0.2)
#'
#' set.seed(123)
#' null_model(d)
#' set.seed(123)
#' null_model(d,null_type="subnetwork1")
#' set.seed(123)
#' null_model(d,null_type="Savue")
#' set.seed(123)
#' null_model(d,number=2,null_type="Savue")


null_model<-function(network, number=NULL, null_type=c("subnetwork1","subnetwork2","all","Savue")){
   mat<-as.matrix(network[])
   mat1<-mat[V(network)$level==0,V(network)$level==1]
   mat2<-mat[V(network)$level==1,V(network)$level==2]
   if(missing(null_type))
      null_type<-"all"
   if(missing(number))
      number<-1
   if(null_type=="subnetwork1"){
      null_list<-vaznullR(number,mat1)
      MAT<-mat
      null_network<-lapply(null_list, function(x){
         MAT[V(network)$level==0,V(network)$level==1]<-x
         Nnetwork<-graph_from_adjacency_matrix(MAT)
         V(Nnetwork)$name<-V(network)$name
         V(Nnetwork)$level<-V(network)$level
         return(Nnetwork)
      })
   }
   else if(null_type=="subnetwork2"){
      null_list<-vaznullR(number,mat2)
      MAT<-mat
      null_network<-lapply(null_list, function(x){
         MAT[V(network)$level==1,V(network)$level==2]<-x
         Nnetwork<-graph_from_adjacency_matrix(MAT)
         V(Nnetwork)$name<-V(network)$name
         V(Nnetwork)$level<-V(network)$level
         return(Nnetwork)
      })
   }
   else if(null_type=="all"){
      null_list1<-vaznullR(number,mat1)
      null_list2<-vaznullR(number,mat2)
      MAT<-mat
      null_network<-lapply(1:number,function(x){
         MAT[V(network)$level==0,V(network)$level==1]<-null_list1[[x]]
         MAT[V(network)$level==1,V(network)$level==2]<-null_list2[[x]]
         Nnetwork<-graph_from_adjacency_matrix(MAT)
         V(Nnetwork)$name<-V(network)$name
         V(Nnetwork)$level<-V(network)$level
         return(Nnetwork)
      })
   }
   else if(null_type=="Savue"){
      null_list1<-SavueR(number,mat1,type="col")
      null_list2<-SavueR(number,mat2,type="row")
      MAT<-mat
      null_network<-lapply(1:number,function(x){
         MAT[V(network)$level==0,V(network)$level==1]<-null_list1[[x]]
         MAT[V(network)$level==1,V(network)$level==2]<-null_list2[[x]]
         Nnetwork<-graph_from_adjacency_matrix(MAT)
         V(Nnetwork)$name<-V(network)$name
         V(Nnetwork)$level<-V(network)$level
         return(Nnetwork)
      })
   }
   else
      stop("Error: null_type is not a valid input!")

   return(null_network)
}
