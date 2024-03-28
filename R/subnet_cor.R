#' Correlation of the Structural indices between Subnetworks
#'
#' Calculating correlation of interconnecting species generalism among Subnetworks ("general_'~'cor") and correlation between similarities of interconnecing species interaction partners in two subnetworks ("similar_'~'cor").
#'
#' @param network.or.subnet_mat1 Either a 'igraph' multilayer(tripartite) network, or a numeric matrix(or data.frame) representing interactions between two groups of species. The network containS
#'  interlayer links and without intralayer links. Each row and column of matrix represents the species in the first and second layers of the tripartite network,respectively. Elements of matrix are non-zero
#'  nnumebers if the interlayer species are connected, and 0 otherwise.
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents the species in the second and third layers of the
#'  tripartite network,respectively. Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise.
#' @param weighted Logical. should elements of matrix be fractional? Default to FALSE. Generally, 'igraph' network represent a spare matrix, so \code{weighted} is FALSE. While elements of matrix represent interaction strength, \code{weighted} is TRUE.
#' @details
#'
#' \strong{weighted}
#'
#' If the \code{weighted} = FALSE, the input for the parameter can be:
#' \code{network.or.subnet_mat1}: input a 'igraph' of network data independently or input sparse matrix together with \code{subnet_mat2}.
#'
#'
#' If the \code{weighted} = TRUE, the input for the parameter can be:
#' \code{network.or.subnet_mat1}: must input matrix(or data.frame) together with \code{subnet_mat2}. the matrix can be sparse matrix and matrix of interaction strength.
#'
#' @return
#' Return a numeric vector of two elements representing correlation of interconnecting species generalism among Subnetworks and correlation between similarities of interconnecing species interaction partners in two subnetworks.
#'
#' If \code{weighted} = FALSE, the results will show "general_cor=  ;similar_cor=  ;" and If \code{weighted} = TRUE, the results will show "general_weight_cor=  ;similar_weight_cor=  ;"
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
#' subnet_cor(d)
#'
#' md1<-matrix(sample(c(0,1),100,replace=TRUE),10,10)
#' md2<-matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' subnet_cor(md1,md2)
#' subnet_cor(md1,md2,weighted=TRUE)
#'
#' mdw1<-matrix(runif(100,0,1),10,10)
#' mdw2<-matrix(runif(120,0,1),10,12)
#' subnet_cor(mdw1,mdw2,weighted=TRUE)
#'
subnet_cor<-function(network.or.subnet_mat1,subnet_mat2=NULL,weighted=FALSE){
   if(!weighted){
      if(inherits(network.or.subnet_mat1,"igraph")==T){
         network<-adject_net(network.or.subnet_mat1)
         mat<-as.matrix(network[])
         mat1<-t(mat[V(network)$level==0,V(network)$level==1])
         mat2<-mat[V(network)$level==1,V(network)$level==2]
      }
      else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T){
         if(nrow(network.or.subnet_mat1)==nrow(subnet_mat2)){
            mat1<-as.matrix(network.or.subnet_mat1)
            mat2<-as.matrix(subnet_mat2)
            logi<-(as.numeric(rowSums(mat1))*as.numeric(rowSums(mat2)))!=0
            mat1<-mat1[logi,]
            mat2<-mat2[logi,]
         }
         else
            stop("please check the row numeber of this 'two matrix'")
      }
      else{
         stop("please check the type of 'network.or.subnet_mat1'")
      }
      general_cor<-Kendall_cor(as.numeric(rowSums(mat1)),as.numeric(rowSums(mat2)))
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
      message(paste0(c("general_cor","similar_cor"),"= ",seq=c(general_cor,similar_cor),";"),"\n")
      return(c(general_cor,similar_cor))
   }
   else{
      if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))){
         if(nrow(network.or.subnet_mat1)==nrow(subnet_mat2)){
            subnet_mat1<-as.matrix(network.or.subnet_mat1)
            subnet_mat2<-as.matrix(subnet_mat2)
            logi<-(as.numeric(rowSums(subnet_mat1))*as.numeric(rowSums(subnet_mat2)))!=0
            subnet_mat1<-subnet_mat1[logi,]
            subnet_mat2<-subnet_mat2[logi,]

            general_weight1<-apply(subnet_mat1,1,function(x){
               if(sum(x)==0){return(0)}
               else{x<-x[x!=0];return(-sum((x/sum(x))*(log(x/sum(x)))))}
            })
            general_weight2<-apply(subnet_mat2,1,function(x){
               if(sum(x)==0){return(0)}
               else{x<-x[x!=0];return(-sum((x/sum(x))*(log(x/sum(x)))))}
            })
            general_weight_cor<-Kendall_cor(general_weight1,general_weight2)
            jaccard_weight1<-NULL
            jaccard_weight2<-NULL
            for(i in 1:(nrow(subnet_mat1)-1)){
               for(j in (i+1):nrow(subnet_mat2)){
                  jaccard_weight1<-c(jaccard_weight1,sum(apply(subnet_mat1[c(i,j),],2,min))/sum(apply(subnet_mat1[c(i,j),],2,max)))
                  jaccard_weight2<-c(jaccard_weight2,sum(apply(subnet_mat2[c(i,j),],2,min))/sum(apply(subnet_mat2[c(i,j),],2,max)))
               }
            }
            similar_weight_cor<-Kendall_cor(jaccard_weight1,jaccard_weight2)
            message(paste0(c("general_weight_cor","similar_weight_cor"),"= ",seq=c(general_weight_cor,similar_weight_cor),";"),"\n")
            return(c(general_weight_cor,similar_weight_cor))
         }
         else{
            stop("please check the row numeber of this 'two matrix'")
         }
      }
      else{
         stop("please check the type of 'network.or.subnet_mat1'")
      }
   }
}
