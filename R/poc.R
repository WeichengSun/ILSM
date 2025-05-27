#' Proportion of connector nodes
#'
#' Calculating the proportion of connector nodes in the shared set of nodes.
#'
#' @param network.or.subnet_mat1 An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.
#' @param subnet_mat2 A matrix representing one subnetwork.
#'
#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#'
#' Two types of inputs \code{network.or.subnet_mat1} can be processed:
#' \itemize{
#' \item An "igraph" object with node attribute 'level' (0 for a-nodes, 1 for b-nodes,2 for c-nodes).
#' \item Or a matrix representing subnetwork P, and must be input with \code{subnet_mat2} representing subnetwork Q.
#' }
#'
#' If the inputs are two matrices, please make sure the rows of
#'  \code{network.or.subnet_mat1} and \code{subnet_mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
#'  names to produce connector nodes. Otherwise, row numbers are assigned to row names and matched. Within the two matrices (P and Q), columns represents a-group nodes and c-group nodes respectively.
#'  Elements in matrices are non-zero values if two nodes are linked with or without weights, and 0 otherwise.
#'
#'
#' @return
#' Return a value representing the proportion of connector species in shared set of species.
#'
#' @import igraph
#' @export
#' @references
#'
#' Battiston, F., Nicosia, V. & Latora, V. (2014) Structural measures for multiplex networks. Physical Review E, 89, 032804.
#'
#' Domínguez-García, V. and Kéfi, S. (2024). The structure and robustness of ecological networks with two interaction types. PLOS Computational Biology, 20(1), e1011770.
#'
#' Guimera, R. & Amaral, L.A.N. (2005) Cartography of complex networks: modules and universal roles. Journal of Statistical Mechanics: Theory and Experiment, 2005, P02001.
#'
#'
#' @examples
#'
#' ## generate a random binary tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#' poc(Net)
#'
#'## empirical network
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
      #message(paste(c("POC"),"=",seq=c(C)),"\n")
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
            stop("Make sure matrices either have no row names or have full row names. No NA")
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
         #message(paste(c("POC"),"=",seq=c(C)),"\n")
         return(C)
      }
      else
         stop("please check the type of 'subnet_mat2' or the row number of this matrix")
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")
}
