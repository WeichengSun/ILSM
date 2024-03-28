#' Counting the degree hub of multilayer network
#'
#' This function counts degree hub that the proportion of interconnecting species serving as the core node of the network degree.
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
#'
#'
#'
#' @return
#' Print a "Hc= ;" and Return a numberic value representing the degree hub of network.
#'
#' @importFrom igraph V
#' @importFrom stats quantile
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
#' Hc(d)
#'
#' md1<-matrix(sample(c(0,1),100,replace=TRUE),10,10)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' Hc(md1,md2)
#'
#' mdw1<-matrix(sample(c(rep(0,40),runif(60,0,1))),10,10)
#' mdw2<-matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' Hc(mdw1,mdw2)
#'
Hc<-function(network.or.subnet_mat1,subnet_mat2=NULL){
   if(inherits(network.or.subnet_mat1,"igraph")){
      network<-network.or.subnet_mat1
      mat<-as.matrix(network[])
      mat1<-t(mat[V(network)$level==0,V(network)$level==1])
      mat2<-mat[V(network)$level==1,V(network)$level==2]
      logi<-rowSums(mat1)*rowSums(mat2)!=0
      m1<-mat1[logi,]
      m2<-mat2[logi,]
      link_degree<-sort(c((rowSums(mat1)+rowSums(mat2)),colSums(mat1),colSums(mat2)),decreasing = T)
      five_summ<-quantile(link_degree,0.8)
      H_C<-sum((rowSums(m1)+rowSums(m2))>=five_summ)/sum(logi)
      message(paste(c("H_C"),"=",seq=c(H_C)),"\n")
      return(H_C)
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))){
      if(inherits(subnet_mat2,c("matrix","data.frame"))&(nrow(network.or.subnet_mat1)==nrow(subnet_mat2))){
         mat1<-as.matrix(network.or.subnet_mat1)
         mat1[mat1>0]<-1
         mat2<-as.matrix(subnet_mat2)
         mat2[mat2>0]<-1
         logi<-rowSums(mat1)*rowSums(mat2)!=0
         m1<-mat1[logi,]
         m2<-mat2[logi,]
         link_degree<-sort(c((rowSums(mat1)+rowSums(mat2)),colSums(mat1),colSums(mat2)),decreasing = T)
         five_summ<-quantile(link_degree,0.8)
         H_C<-sum((rowSums(m1)+rowSums(m2))>=five_summ)/sum(logi)
         message(paste(c("H_C"),"=",seq=c(H_C)),"\n")
         return(H_C)
      }
      else{
         stop("please check the type of 'subnet_mat2' or the row numeber of this 'two matrix'")
      }
   }
   else{
      stop("please check the type of 'network.or.subnet_mat1'")
   }
}
quantile
