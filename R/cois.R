#' Correlation of Interaction Similarity for Shared species: CoIS
#'
#' Calculating correlation of interaction similarity for shared species("COis_ ~") in two subnetworks.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class which contains three groups of species and interactions within layers without interactions between each group of species, or a numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#'
#' @param weighted Logical. should elements of matrix be fractional? Default to FALSE. Generally, 'igraph' network represent a spare matrix, so \code{weighted} is FALSE. While elements of matrix represent interaction strength, \code{weighted} is TRUE.
#' @details
#'
#' \strong{weighted}
#'
#' If the \code{weighted} = FALSE, the input for the parameter can be:
#' \itemize{\item{\code{network.or.subnet_mat1}: input a 'igraph' of network data independently or input sparse matrix together with \code{subnet_mat2}.}}
#'
#' If the \code{weighted} = TRUE, the input for the parameter can be:
#' \itemize{\item{\code{network.or.subnet_mat1}: must input matrix(or data.frame) together with \code{subnet_mat2}. the matrix can be sparse matrix and matrix of interaction strength.}}
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
#'
#' \item{When the two matrices can have different numbers of rows:}
#' \itemize{
#' \item{(1). If both matrices have row names, then the function counts all row names to produce two new matrices with the same row names.}
#' \item{(2). If at most one matrix has row names, the function assigns new row names to both matrices on a row-to-row basis (any extra row names are assigned a new value) and then counts all row names to produce two new matrices with the same row names.}
#' }
#'
#' \item{When the two matrices can have the same numbers of rows:}
#' \itemize{
#' \item{No matter how the row names of the two matrices are arranged, as long as the row names are exactly the same; But we don't handle matrices with empty row names (the function will give an error).}
#' }
#'
#' \item{The two matrices can have different numbers of rows, but read our default handling carefully to make sure the calculation is accurate when using this function!!!}
#' }
#' About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}.
#'
#'
#' @return
#' Return a numeric value representing correlation of interaction similarity for shared species among subnetworks .
#'
#' If \code{weighted} = FALSE, the results will show "CoIS= ;" and If \code{weighted} = TRUE, the results will show "CoIS_weight= ;"
#' @import igraph
#' @export
#' @references
#'
#' Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
#'
#'
#' @examples
#'
#' set.seed(12)
#' d <- build_net(11,15,17,0.2)
#' cois(d)
#'
#' md1<-matrix(sample(c(0,1),110,replace=TRUE),10,11)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' cois(md1,md2)
#' cois(md1,md2,weighted=TRUE)
#'
#' md1<-matrix(sample(c(0,1),80,replace=TRUE),8,10)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' cois(md1,md2)
#'
#' mdw1<-matrix(runif(110,0,1),10,11)
#' mdw2<-matrix(runif(120,0,1),10,12)
#' cois(mdw1,mdw2,weighted=TRUE)
#'
#' set.seed(1)
#' mdw1<-matrix(runif(80,0,1),8,10)
#' mdw2<-matrix(runif(120,0,1),10,12)
#' cois(mdw1,mdw2,weighted=TRUE)
#'
#'

cois<-function(network.or.subnet_mat1, subnet_mat2=NULL, weighted=FALSE){
   if(!weighted){
      if(inherits(network.or.subnet_mat1,"igraph")==T){
         network<-adject_net(network.or.subnet_mat1)
         mat<-as.matrix(network[])
         mat1<-t(mat[V(network)$level==0,V(network)$level==1])
         mat2<-mat[V(network)$level==1,V(network)$level==2]
      }
      else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
         mat1<-network.or.subnet_mat1
         mat2<-subnet_mat2
         if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
            rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
            rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         }
         if(nrow(mat1)!=nrow(mat2))
            message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
         if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         else
            stop("Make sure matrices either have no row names or have full row names. No NA!!!")
         mat_1<-matrix(0,length(matrow),ncol(mat1))
         rownames(mat_1)<-matrow
         mat_1[rownames(mat1),]<-mat1
         mat_1[mat_1>0]<-1
         mat_2<-matrix(0,length(matrow),ncol(mat2))
         rownames(mat_2)<-matrow
         mat_2[rownames(mat2),]<-mat2
         mat_2[mat_2>0]<-1
         mat1<-mat_1
         mat2<-mat_2
      }
      else
         stop("please check the type of 'network.or.subnet_mat1'")
      logi<-(as.numeric(rowSums(mat1))*as.numeric(rowSums(mat2)))!=0
      mat1<-mat1[logi,]
      mat2<-mat2[logi,]
      jaccard_vector1<-NULL
      jaccard_vector2<-NULL
      for(i in 1:(nrow(mat1)-1)){
         for(j in (i+1):nrow(mat2)){
            same_degree1<-sum(mat1[i,]*mat1[j,])
            sum_degree1<-sum((mat1[i,]+mat1[j,])>0)
            jaccard_vector1<-c(jaccard_vector1,same_degree1/sum_degree1)
            same_degree2<-sum(mat2[i,]*mat2[j,])
            sum_degree2<-sum((mat2[i,]+mat2[j,])>0)
            jaccard_vector2<-c(jaccard_vector2,same_degree2/sum_degree2)
         }
      }
      jaccard_vector1[is.na(jaccard_vector1)]<-0
      jaccard_vector2[is.na(jaccard_vector2)]<-0
      similar_cor<-Kendall_cor(jaccard_vector1,jaccard_vector2)
      message(paste0("CoIS= ",seq=similar_cor,";"),"\n")
      return(similar_cor)
   }
   else{
      if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
         mat1<-network.or.subnet_mat1
         mat2<-subnet_mat2
         if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
            rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
            rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         }
         if(nrow(mat1)!=nrow(mat2))
            message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
         if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         else
            stop("Make sure matrices either have no row names or have full row names. No NA!!!")
         mat_1<-matrix(0,length(matrow),ncol(mat1))
         rownames(mat_1)<-matrow
         mat_1[rownames(mat1),]<-mat1
         mat_2<-matrix(0,length(matrow),ncol(mat2))
         rownames(mat_2)<-matrow
         mat_2[rownames(mat2),]<-mat2
         mat1<-mat_1
         mat2<-mat_2
      }
      else
         stop("please check the type of 'network.or.subnet_mat1'")
      subnet_mat1<-mat1
      subnet_mat2<-mat2
      logi<-(as.numeric(rowSums(subnet_mat1))*as.numeric(rowSums(subnet_mat2)))!=0
      subnet_mat1<-subnet_mat1[logi,]
      subnet_mat2<-subnet_mat2[logi,]

      jaccard_weight1<-NULL
      jaccard_weight2<-NULL
      for(i in 1:(nrow(subnet_mat1)-1)){
         for(j in (i+1):nrow(subnet_mat2)){
            jaccard_weight1<-c(jaccard_weight1,sum(apply(subnet_mat1[c(i,j),],2,min))/sum(apply(subnet_mat1[c(i,j),],2,max)))
            jaccard_weight2<-c(jaccard_weight2,sum(apply(subnet_mat2[c(i,j),],2,min))/sum(apply(subnet_mat2[c(i,j),],2,max)))
         }
      }
      similar_weight_cor<-Kendall_cor(jaccard_weight1,jaccard_weight2)
      message(paste0("CoIS_weight= ",seq=similar_weight_cor,";"),"\n")
      return(similar_weight_cor)
   }
}
