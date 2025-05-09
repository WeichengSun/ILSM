#' Correlation of degree of connector nodes: CoID
#'
#' Calculating correlation of degree of connector nodes ("CoID_~").
#'
#' @param network.or.subnet_mat1 An igraph object or matrix. A tripartite network contains three groups (a,b,c) of nodes and two subnetworks (P including a- and b-groups of species, Q including b- and c-groups of species).If the network is weighted, "matrix" is recommended, otherwise unweighted result is returned. If matrices are provided, this should be the matrix representing subnetwork P that have links between a- and b-group species.
#'  Within this matrix, each row and column represents b-group species and a-group species respectively. For matrix, row names are required, otherwise row number is used.
#'  Elements of matrix are non-zero values if two species are linked, and 0 otherwise.
#'
#' @param subnet_mat2 The matrix representing subnetwork Q that have links between b- and c-group species.
#'  Within this matrix, each row and column represents b-group species and c-group species respectively.
#'  Elements of matrix are non-zero values if two species are linked, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#'
#' @param weighted Logical. If true, shannon diveristy of interaction strengths of connector nodes are used. Default to FALSE.
#' @param weight_type For weighted networks,the definition of weighted degree for a connector node, supporting "shannon" or "sum"
#' @param method  Correlation method ("pearson", "kendall" or "spearman"). Default to "kendall".
#' @details
#'This function follows Sauve et al.(2016)) to calculate the correlation of interaction degree (or Shannon diversity of interaction strength ) of connector nodes. For the binary network, connector nodes' degree is calculated in each subnetwork.
#'For the quantitative network, Shannon diversity of interaction strength for each connector node (i) is calculated as ??. V is the set of a or c nodes it interacts.
#'Three correlation methods are provided. Kendall correlation is recommended following Sauve et al.(2016).
#' \strong{weighted}
#'
#' If the \code{weighted} = FALSE, the input network can be an "igraph" object or two matrices.If a weighted network is provided, it will be transformed to a binary network.
#' If the \code{weighted} = TRUE, the input network can only be two matrices. Correlation of Shannon diversity of interaction strengths for connector nodes are returned.
#'
#' \strong{network.or.subnet_mat1 and subnet_mat2}
#'
#' There are two types of inputs \code{network.or.subnet_mat1} that can be processed:
#' \itemize{
#' \item Input is a network of type "igraph" alone.
#' \item Must be entered as matrix with \code{subnet_mat2}.
#' }
#'
#' If the inputs are two matrices, please make sure the rows of
#'  \code{network.or.subnet_mat1} and \code{subnet_mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
#' names to produce connector nodes. Otherwise, row numbers are assigned as row names.
#'
#' @return
#' Return a numeric value representing correlation of interaction degree: CoID.
#'
#' If \code{weighted} = FALSE, the results will show "CoID= ;" and If \code{weighted} = TRUE, the results will show "CoID_weight= ;"
#' @import igraph
#' @export
#'
#' @srrstats {G2.14b,G2.16} It ignores to handle missing (`NA`) data.
#'
#' @references
#'
#' Sauve, A. M., ThC)bault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
#'
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#' coid(Net)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' coid(Net)
#'
#'##input as binary matrices,with row names.
#' md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
#' md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
#' coid(md1,md2)
#' coid(md1,md2,weighted=T)
#'
#'##input as weighted matrices,with row numbers as row names.
#' mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
#' mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#'coid(mdw1,mdw2)
#'coid(mdw1,mdw2,weighted=T)
#'


coid<-function(network.or.subnet_mat1, subnet_mat2=NULL, weighted=FALSE,weight_type="shannon",method="kendall" ){
   if(!weighted){
      if(inherits(network.or.subnet_mat1,"igraph")==T){
         network<-adjust_net(network.or.subnet_mat1)
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
         stop("please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
      logi<-(as.numeric(rowSums(mat1))*as.numeric(rowSums(mat2)))!=0
      mat1<-mat1[logi,]
      mat2<-mat2[logi,]
      general_cor<-cor(as.numeric(rowSums(mat1)),as.numeric(rowSums(mat2)), method=method )
      message(paste0("CoID= ",seq=round(general_cor,8),";"),"\n")
      return(general_cor)
   }
   else{
      if(inherits(network.or.subnet_mat1,c("matrix"))==T && inherits(subnet_mat2,c("matrix"))==T){
         mat1<-network.or.subnet_mat1
         mat2<-subnet_mat2
         if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
            rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
            rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         }
         if(nrow(mat1)!=nrow(mat2))
            message("Different row numbers for the two matrices, and connector nodes are determined by matching row names.")
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
         stop("please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
      subnet_mat1<-mat1
      subnet_mat2<-mat2
      logi<-(as.numeric(rowSums(subnet_mat1))*as.numeric(rowSums(subnet_mat2)))!=0
      subnet_mat1<-subnet_mat1[logi,]
      subnet_mat2<-subnet_mat2[logi,]

      general_weight1<-apply(subnet_mat1,1,function(x){
         if(sum(x)==0){return(0)}
         else{x<-x[x!=0];
         if (weight_type="shannon"){return(-sum((x/sum(x))*(log(x/sum(x))))}
         else if (weight_type="sum"){return(sum(x))}
         else{ stop("weight_type should be 'shannon' or 'sum'")}
         )}
      })
      general_weight2<-apply(subnet_mat2,1,function(x){
         if(sum(x)==0){return(0)}
         else{x<-x[x!=0];
         if (weight_type="shannon"){return(-sum((x/sum(x))*(log(x/sum(x))))}
         else if (weight_type="sum"){return(sum(x))}
         else{ stop("weight_type should be 'shannon' or 'sum'")}
         )}
      })
      general_weight_cor<-cor(general_weight1,general_weight2,method=method)
      message(paste0("CoID_weight= ",seq=round(general_weight_cor,8),";"),"\n")
      return(general_weight_cor)
   }
}
