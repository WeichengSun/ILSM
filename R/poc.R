#' Proportion of connector nodes
#'
#' Calculating the proportion of connector nodes in the shared set of nodes.
#'
#' @param network.or.subnet_mat1 Igraph or matrix. Both weighted or unweighted networks are accepted. A tripartite network of 'igraph' class which contains three groups (a,b,c) of species and two subnetworks (P including a and b groups of species, Q including b and c groups of species). If matrices are provided, this should be the matrix representing subnetwork P that have links between a- and b-group species. 
#'  Within this matrix, each row and column represents b-group species and a-group species respectively. For matrix, row names are required, otherwise row number is used. 
#'  Elements of matrix are non-zero values if two species are linked, and 0 otherwise.
#' @param subnet_mat2 The matrix representing subnetwork Q that have links between b- and c-group species. 
#'  Within this matrix, each row and column represents b-group species and c-group species respectively.
#'  Elements of matrix are non-zero values if two species are linked, and 0 otherwise.If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#' @details
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
#'
#' @return
#' Print a "POC= ;" and Return a numeric value representing the proportion of connector species in shared group of species.
#'
#' @import igraph
#' @export
#' @references
#'
#' Battiston, F., Nicosia, V. & Latora, V. (2014) Structural measures for multiplex networks. Physical Review E, 89, 032804.
#'
#' Domínguez-García, V., & Kéfi, S. (2024). The structure and robustness of ecological networks with two interaction types. PLOS Computational Biology, 20(1), e1011770.
#'
#' Guimera, R. & Amaral, L.A.N. (2005) Cartography of complex networks: modules and universal roles. Journal of Statistical Mechanics: Theory and Experiment, 2005, P02001.
#'
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' poc(Net)
#'
#'##input as binary matrices,with row names.
#' md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
#' md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
#' poc(md1,md2)
#' 
#'##input as weighted matrices,with row numbers as row names.
#' mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
#' mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' poc(mdw1,mdw2)
#'
#'
poc<-function(network.or.subnet_mat1, subnet_mat2=NULL){
   if(inherits(network.or.subnet_mat1,"igraph")){
      network<-network.or.subnet_mat1
      mat<-as.matrix(network[])
      mat1<-t(mat[V(network)$level==0,V(network)$level==1])
      mat2<-mat[V(network)$level==1,V(network)$level==2]
      logi<-rowSums(mat1)*rowSums(mat2)!=0
      C<-sum(logi)/nrow(mat1)
      message(paste(c("POC"),"=",seq=c(C)),"\n")
      return(C)
   }
   else if(inherits(network.or.subnet_mat1,c("matrix"))){
      if(inherits(subnet_mat2,c("matrix"))){
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
         C<-sum(logi)/nrow(mat1)
         message(paste(c("POC"),"=",seq=c(C)),"\n")
         return(C)
      }
      else
         stop("please check the type of 'subnet_mat2' or the row number of this matrix")
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")
}
