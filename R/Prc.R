#' Participation ratio of interconnecting code
#'
#' Counting participation ratio that the difference in the degree of interconnecting nodes within two subnetworks.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class which contains interlayer links and without intralayer links, or a numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if two groups of species are connected, and 0 otherwise.
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
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
#' Print a "PR_C= ;" and Return a numeric value representing difference in the degree of interconnecting nodes within two subnetworks of multialyer network.
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
#' Prc(d)
#'
#' md1<-matrix(sample(c(0,1),80,replace=TRUE),8,10)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' Prc(md1,md2)
#'
#' mdw1<-matrix(sample(c(rep(0,60),runif(60,0,1))),12,10)
#' mdw2<-matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' Prc(mdw1,mdw2)
#'
Prc<-function(network.or.subnet_mat1, subnet_mat2=NULL){
   if(inherits(network.or.subnet_mat1,"igraph")){
      network<-network.or.subnet_mat1
      mat<-as.matrix(network[])
      mat1<-t(mat[V(network)$level==0,V(network)$level==1])
      mat2<-mat[V(network)$level==1,V(network)$level==2]
      logi<-rowSums(mat1)*rowSums(mat2)!=0
      m1<-mat1[logi,]
      m2<-mat2[logi,]
      PR_C<-mean(2*apply(rbind(rowSums(m1),rowSums(m2)),2,min)/apply(rbind(rowSums(m1),rowSums(m2)),2,sum))
      message(paste(c("PR_C"),"=",seq=c(PR_C)),"\n")
      return(PR_C)
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
         m1<-mat1[logi,]
         m2<-mat2[logi,]
         PR_C<-mean(2*apply(rbind(rowSums(m1),rowSums(m2)),2,min)/apply(rbind(rowSums(m1),rowSums(m2)),2,sum))
         message(paste(c("PR_C"),"=",seq=c(PR_C)),"\n")
         return(PR_C)
      }
      else
         stop("please check the type of 'subnet_mat2' or the row numeber of this 'two matrix'")
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")
}
