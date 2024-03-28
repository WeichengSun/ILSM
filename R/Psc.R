#' Proportion of interconnection
#'
#' Calculating the proportion of species sharing with other species of two subnetworks in intermediate layer.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class, or a numeric matrix(or data.frame) representing interactions between two groups of species. The network contains
#'  interlayer links and without intralayer links. Each row and column of matrix represents the species in the first and second layers of the tripartite network respectively.
#'  Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise.
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents the species in the second and third layers of the
#'  tripartite network respectively. Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#' @details
#'
#' \strong{network.or.subnet_mat1}
#' There are two types of data that can be processed:
#' (1). Input in a network of type "igraph" alone
#' (2). Must be entered as data frame or matrix with \code{subnet_mat2}
#'
#' About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}
#'
#' @return
#' Print a "PS_C= ;" and Return a numberic value representing the proportion of sharing species in intermediate layer.
#'
#' @import igraph
#' @export
#' @references
#'
#' Domínguez-García, V., & Kéfi, S. (2021). The structure and robustness of tripartite ecological networks. bioRxiv, 2021-10.
#'
#'
#' @examples
#'
#' set.seed(15)
#' d <- build_net(11,15,17,0.2)
#' Psc(d)
#'
#' md1<-matrix(sample(c(0,1),100,replace=TRUE),10,10)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' Psc(md1,md2)
#'
#' mdw1<-matrix(sample(c(rep(0,40),runif(60,0,1))),10,10)
#' mdw2<-matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' Psc(mdw1,mdw2)
#'
Psc<-function(network.or.subnet_mat1,subnet_mat2=NULL){
   if(inherits(network.or.subnet_mat1,"igraph")){
      network<-network.or.subnet_mat1
      mat<-as.matrix(network[])
      mat1<-t(mat[V(network)$level==0,V(network)$level==1])
      mat2<-mat[V(network)$level==1,V(network)$level==2]
      logi<-rowSums(mat1)*rowSums(mat2)!=0
      C<-sum(logi)/nrow(mat1)
      message(paste(c("C"),"=",seq=c(C)),"\n")
      return(C)
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))){
      if(inherits(subnet_mat2,c("matrix","data.frame"))&(nrow(network.or.subnet_mat1)==nrow(subnet_mat2))){
         mat1<-as.matrix(network.or.subnet_mat1)
         mat1[mat1>0]<-1
         mat2<-as.matrix(subnet_mat2)
         mat2[mat2>0]<-1
         logi<-rowSums(mat1)*rowSums(mat2)!=0
         C<-sum(logi)/nrow(mat1)
         message(paste(c("C"),"=",seq=c(C)),"\n")
         return(C)
      }
      else{
         stop("please check the type of 'subnet_mat2' or the row numeber of this 'two matrix'")
      }
   }
   else{
      stop("please check the type of 'network.or.subnet_mat1'")
   }
}
