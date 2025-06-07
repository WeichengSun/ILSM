#' Correlation of Interaction Similarity for connector nodes
#'
#' Calculating correlation of interaction similarity for connector nodes.
#'
#' @param network.or.subnet_mat1 An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.
#' @param subnet_mat2 A matrix representing one subnetwork.
#' @param weighted Logical. Default to FALSE. See details.
#' @param method  Correlation method ("pearson", "kendall" or "spearman"). Default to "kendall".
#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#'
#' This function follows Sauve et al.(2016) to calculate the correlation of interaction similarity of connector nodes. When weighted=FALSE, the Jaccard similarity coefficients between all pairs of connector nodes in each subnetwork, defined as the number of shared interaction partners divided
#' by the total number of interaction partners involved (i.e. the maximum number of interactions they could share). When weighted=TRUE, the quantitative
#' similarities of all pairs of connector nodes in each subnetwork were measured using the generalised Jaccard coefficient.
#' Three correlation methods are provided. Kendall correlation is recommended following Sauve et al.(2016).
#'
#' Two types of inputs \code{network.or.subnet_mat1} can be processed:
#' \itemize{
#' \item An "igraph" object with node attribute 'level' (0 for a-nodes, 1 for b-nodes,2 for c-nodes). If the input is a weighted network, the edge should have a 'weight' attribute.
#' \item Or a matrix representing subnetwork P, and must be input with \code{subnet_mat2} representing subnetwork Q.
#' }
#'
#' If the inputs are two matrices, please make sure the rows of
#'  \code{network.or.subnet_mat1} and \code{subnet_mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
#'  names to produce connector nodes. Otherwise, row numbers are assigned to row names and matched. Within the two matrices (P and Q), columns represents a-group nodes and c-group nodes respectively.
#'  Elements in matrices are non-zero values if two nodes are linked with or without weights, and 0 otherwise.
#'
#' @return
#' Return a numeric value representing correlation of interaction similarity for connector species.
#' @import igraph
#' @export
#' @references
#'
#' Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
#'
#'
#' @srrstats {G1.5} This standard belongs here.
#' @srrstats {G2.2} This standard belongs here.
#'
#'
#' @examples
#'
#' ## generate a random binary tripartite network
#' set.seed(12)
#' Net <- build_toy_net(11,15,16,0.2)
#' cois(Net)

#' ## empirical network
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' cois(Net)
#' set.seed(13)
#' library(igraph)
#' E(Net)$weight<-runif(length(E(Net)),0.1,1.5)#random weights assigned
#' cois(Net,weighted=T)
#'
#'##input as binary matrices,with row names.
#' md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
#' md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
#' cois(md1,md2)
#'
#'##input as weighted matrices,with row numbers as row names.
#' mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
#' mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#'cois(mdw1,mdw2)
#'cois(mdw1,mdw2,weighted=T)
#'

cois<-function(network.or.subnet_mat1, subnet_mat2=NULL, weighted=FALSE,method="kendall"){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network<-adjust_net(network.or.subnet_mat1,weighted=T)
      mat<-as.matrix(network[])
      mat1<-t(mat[V(network)$level==0,V(network)$level==1])
      mat2<-mat[V(network)$level==1,V(network)$level==2]
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
      mat1<-network.or.subnet_mat1
      mat2<-subnet_mat2
      if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
         message("No rownames for matrices, so row IDs are used!")
         rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
         rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
         matrow<-unique(c(rownames(mat1),rownames(mat2)))
      }
      #if(nrow(mat1)!=nrow(mat2))
      #   message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
      if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0){
         matrow<-unique(c(rownames(mat1),rownames(mat2)))
         if(length(matrow)==0) stop("No connectors existed.")
      }
      else {stop("Please make sure matrices either have no row names or have full row names. No NA!")}

      mat_1<-matrix(0,length(matrow),ncol(mat1))
      rownames(mat_1)<-matrow
      mat_1[rownames(mat1),]<-mat1
      #mat_1[mat_1>0]<-1
      mat_2<-matrix(0,length(matrow),ncol(mat2))
      rownames(mat_2)<-matrow
      mat_2[rownames(mat2),]<-mat2
      #mat_2[mat_2>0]<-1
      mat1<-mat_1
      mat2<-mat_2
   }
   else
      stop("Please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")

   if(!weighted){
      mat1[mat1>0]<-1
      mat2[mat2>0]<-1

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
      similar_cor<-cor(jaccard_vector1,jaccard_vector2,method=method)
      #message(paste0("CoIS= ",seq=similar_cor,";"),"\n")
      return(similar_cor)
   } else{

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
      similar_weight_cor<-cor(jaccard_weight1,jaccard_weight2,method=method)
      #message(paste0("CoIS_weighted= ",seq=similar_weight_cor,";"),"\n")
      return(similar_weight_cor)
   }
}
