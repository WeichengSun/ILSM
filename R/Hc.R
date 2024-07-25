#' Counting the degree hub of multilayer network
#'
#' This function counts degree hub that the proportion of interconnecting species serving as the core node of the network degree.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class which contains three groups of species and interactions within layers without interactions between each group of species, or a numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#'
#' @details
#'
#' \strong{network.or.subnet_mat1} and \strong{subnet_mat2}
#'
#' There are two types of \code{network.or.subnet_mat1} that can be processed:
#' \itemize{
#' \item{(1). Input in a network of type "igraph" alone.}
#' \item{(2). Must be entered as data frame or matrix with \code{subnet_mat2}.}
#' }
#'
#' If the type of inputting is data frame or matrix, please make sure the row of \code{network.or.subnet_mat1} and \code{subnet_mat2} correspond with the second group of species that both belong to two subnetworks and interact with other groups of species.
#' \itemize{
#' \item{Try to make the rows of both matrices have the same attributes. Or we default:}
#' \item{(1). If both matrices have row names, then the function counts all row names to produce two new matrices with the same row names.}
#' \item{(2). If at most one matrix has row names, the function assigns new row names to both matrices on a row-to-row basis (any extra row names are assigned a new value) and then counts all row names to produce two new matrices with the same row names.}
#' \item{The two matrices can have different numbers of rows, but read our default handling carefully to make sure the calculation is accurate when using this function!!!}
#' }
#'
#' About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}.
#'
#'
#' @return
#' Print a "hc= ;" and Return a numeric value representing the degree hub of network.
#'
#' @importFrom igraph V
#' @export
#' @references
#'
#' Battiston, F., Nicosia, V. & Latora, V. (2014) Structural measures for multiplex networks. Physical Review E, 89, 032804.
#'
#' Domínguez-García, V., & Kéfi, S. (2024). The structure and robustness of ecological networks with two interaction types. PLOS Computational Biology, 20(1), e1011770.
#'
#' Guimera, R. & Amaral, L.A.N. (2005) Cartography of complex networks: modules and universal roles. Journal of Statistical Mechanics: Theory and Experiment, 2005, P02001.
#'
#' @examples
#'
#' set.seed(15)
#' d <- build_net(11,15,17,0.2)
#' hc(d)
#'
#' md1<-matrix(sample(c(0,1),80,replace=TRUE),8,10)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' hc(md1,md2)
#'
#' mdw1<-matrix(sample(c(rep(0,60),runif(60,0,1))),12,10)
#' mdw2<-matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' hc(mdw1,mdw2)
#'
hc<-function(network.or.subnet_mat1, subnet_mat2=NULL){
   if(inherits(network.or.subnet_mat1,"igraph")){
      network<-network.or.subnet_mat1
      mat<-as.matrix(network[])
      mat1<-t(mat[V(network)$level==0,V(network)$level==1])
      mat2<-mat[V(network)$level==1,V(network)$level==2]
      logi<-rowSums(mat1)*rowSums(mat2)!=0
      m1<-mat1[logi,]
      m2<-mat2[logi,]
      # link_degree<-sort(c((rowSums(mat1)+rowSums(mat2)),colSums(mat1),colSums(mat2)),decreasing = T)
      link_degree<-order(rowSums(mat1)+rowSums(mat2),decreasing = T)
      hub<-sum(logi[link_degree[1:round(length(link_degree)*0.2)]])
      # five_summ<-quantile(link_degree,0.8)
      # H_C<-sum((rowSums(m1)+rowSums(m2))>=five_summ)/sum(logi)
      H_C<-hub/round(length(link_degree)*0.2)
      message(paste(c("hc"),"=",seq=c(H_C)),"\n")
      return(H_C)
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))){
      if(inherits(subnet_mat2,c("matrix","data.frame"))){
         if(is.null(rownames(network.or.subnet_mat1)) | is.null(rownames(subnet_mat2))){
            rownames(network.or.subnet_mat1)<-paste0("spe",seq=1:nrow(network.or.subnet_mat1))
            rownames(subnet_mat2)<-paste0("spe",seq=1:nrow(subnet_mat2))
            matrow<-unique(c(rownames(network.or.subnet_mat1),rownames(subnet_mat2)))
         }
         if(!is.null(rownames(network.or.subnet_mat1)) & !is.null(rownames(subnet_mat2)) & sum(is.na(rownames(network.or.subnet_mat1)))==0 & sum(is.na(rownames(subnet_mat2)))==0)
            matrow<-unique(c(rownames(network.or.subnet_mat1),rownames(subnet_mat2)))
         else
            stop("Make sure matrices either have no row names or have full row names. No NA!!!")
         mat1<-matrix(0,length(matrow),ncol(network.or.subnet_mat1))
         rownames(mat1)<-matrow
         mat1[rownames(network.or.subnet_mat1),]<-network.or.subnet_mat1
         mat1[mat1>0]<-1
         mat2<-matrix(0,length(matrow),ncol(subnet_mat2))
         rownames(mat2)<-matrow
         mat2[rownames(subnet_mat2),]<-subnet_mat2
         mat2[mat2>0]<-1

         logi<-rowSums(mat1)*rowSums(mat2)!=0
         link_degree<-order(rowSums(mat1)+rowSums(mat2),decreasing = T)
         hub<-sum(logi[link_degree[1:round(length(link_degree)*0.2)]])
         H_C<-hub/round(length(link_degree)*0.2)
         message(paste(c("hc"),"=",seq=c(H_C)),"\n")
         return(H_C)


      }
      else
         stop("please check the type of 'subnet_mat2' or the row numeber of this 'two matrix'")
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")
}


