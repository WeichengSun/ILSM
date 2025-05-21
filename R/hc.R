#' The proportion of hub connectors in the shared sets of nodes between subnetworks
#'
#' This function calculates the proportion of shared node hubs that are connectors (HC), i.e. the top x% of the nodes in the shared set of nodes with the highest degree that are connector nodes.
#'
#' @param network.or.subnet_mat1 Igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork.
#' @param subnet_mat2 The matrix representing one subnetwork.
#' @param hubs The top x% of the nodes in the shared set of nodes with the highest degree. default to 20%.
#' @param weighted Logical. Default to FALSE. If TRUE, weighted degree defined as the sum of interaction strengths attached to a connector node is used.

#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.

#' The function counts the proportion of connector nodes in connector node hubs (HC). The connector node hubs are the top x% of shared nodes with the highest degree (Dominguez-Garcia and Kefi 2024). The default x% is 20%. It always equals 1 if all b-nodes are connector nodes.
#'
#' Two types of inputs \code{network.or.subnet_mat1} can be processed:
#' \itemize{
#' \item An "igraph" object with node attribute 'level' (0 for a-nodes, 1 for b-nodes,2 for c-nodes)
#' \item Or a matrix representing subnetwork P, and must be input with \code{subnet_mat2} representing subnetwork Q.
#' }
#'
#' If the inputs are two matrices, please make sure the rows of
#'  \code{network.or.subnet_mat1} and \code{subnet_mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
#'  names to produce connector nodes. Otherwise, row numbers are assigned to row names and matched. Within the two matrices (P and Q), columns represents a-group nodes and c-group nodes respectively.
#'  Elements in matrices are non-zero values if two nodes are linked with or without weights, and 0 otherwise.
#'
#' @return
#' Return a numeric value representing the proportion of connector nodes in node hubs.
#' network.
#'
#' @importFrom igraph V
#' @export
#' @references
#'
#' Battiston, F., Nicosia, V. & Latora, V. (2014) Structural measures for
#' multiplex networks. Physical Review E, 89, 032804.
#'
#' DomC-nguez-GarcC-a, V., & KC)fi, S. (2024). The structure and robustness of
#' ecological networks with two interaction types. PLOS Computational Biology,
#' 20(1), e1011770.
#'
#' Guimera, R. & Amaral, L.A.N. (2005) Cartography of complex networks: modules
#' and universal roles. Journal of Statistical Mechanics: Theory and
#' Experiment, 2005, P02001.
#'
#' @examples
#'
#' ## generate a random binary tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#' hc(Net)
#'
#' empirical network
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' hc(Net)
#' set.seed(13)
#' E(Net)$weight<-runif(length(E(Net)),0.1,1)#random weights assigned
#' pc(Net,weighted=T)
#'
#'##input as binary matrices,with row names.
#' set.seed(12)
#' md1 <- matrix(sample(c(0,1),5*8,replace=TRUE),5,8,dimnames = list(paste0("b",1:5),paste0("c",1:8)))
#' md2 <- matrix(sample(c(0,1),20*30,replace=TRUE),20,30,dimnames = list(paste0("b",1:20),paste0("a",1:30)))
#' hc(md1,md2)

#'
#'##input as weighted matrices,with row numbers as row names.
#' set.seed(17)
#' mdw1 <- matrix(sample(c(rep(0,20),runif(20,0,1))),5,8)
#' mdw2 <- matrix(sample(c(rep(0,500),runif(100,0,1))),20,30)
#' hc(mdw1,mdw2)
#' hc(mdw1,mdw2,weighted=T)
#'


hc<-function(network.or.subnet_mat1, subnet_mat2 = NULL, hubs = 0.2 ,weighted = FALSE){
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
      else {stop("Make sure matrices either have no row names or have full row names. No NA!")}

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

      logi <- (as.numeric(rowSums(mat1)) * as.numeric(rowSums(mat2))) != 0
      link_degree <- order(rowSums(mat1) + rowSums(mat2), decreasing = T)
      HU <- round(length(link_degree) * hubs)
      hub <- sum(logi[link_degree[1:HU]])
      H_C <- hub / HU
      message(paste(c("HC"), "=", seq = c(H_C)), "\n")
      return(H_C)
   }
   else{
      subnet_mat1 <- mat1
      subnet_mat2 <- mat2
      produt_degree = as.numeric(rowSums(subnet_mat1)) * as.numeric(rowSums(subnet_mat2))
      logi <- (produt_degree) != 0
      link_degree <- order(rowSums(mat1) + rowSums(mat2), decreasing = T)
      HU <- round(length(link_degree) * hubs)
      hub <- sum(logi[link_degree[1:HU]])
      H_C <- hub / HU
      message(paste(c("HC_weighted"), "=", seq = c(H_C)), "\n")
      return(H_C)
   }
}
