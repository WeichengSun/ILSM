#' Transforming two matrices into an igraph object.
#'
#' Transform two matrices into an igraph object.
#'
#' @param mat1 A numeric matrix representing the first subnetwork. Rows should be the shared set of species.
#'
#' @param mat2 A numeric matrix representing the second subnetwork. Rows should be the shared set of species.
#'
#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#' Please make sure the rows of \code{mat1} and \code{mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
#' names to define connector nodes. Otherwise, row numbers are assigned to row names and matched, which might produce an incorrected network. Within the two matrices (P and Q), columns represents a-group nodes and c-group nodes respectively.
#' Elements in matrices are non-zero values if two nodes are linked with or without weights, and 0 otherwise.
#
#' @return
#' Return a network of type "igraph".
#' @export
#'
#' @srrstats {G2.8} It makes the matrices to \code{igraph} data.
#' @srrstats {G2.11,G2.12} The parameters \code{mat1} and \code{mat2} have defined columns.
#' @srrstats {G5.2,G5.2a} Here the message produced within R code by 'stop()' and 'message()'.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph V
#'
#' @examples
#'
#' set.seed(12)
#' MAT <- build_toy_net(11,22,21,0.2,output_matrices=T)
#' M <- trigraph_from_mat(t(MAT[[3]]),MAT[[4]])
#' M

trigraph_from_mat<-function(mat1, mat2, weighted=F){
   if(!inherits(mat1,c("matrix"))|!inherits(mat1,c("matrix"))){
    stop("Please input matrices.")
   }

   if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
    message("Warning! Row IDs were set as rownames for matching connector nodes since no rownames are provided for the matrices")
      rownames(mat1)<-paste0("b",seq=1:nrow(mat1))
      rownames(mat2)<-paste0("b",seq=1:nrow(mat2))
      matrow<-unique(c(rownames(mat1),rownames(mat2)))
   }
   if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
       {matrow<-unique(c(rownames(mat1),rownames(mat2)))}else
      {stop("Please make sure the two matrices have appropriate row names. NA is not accepted.")}

   if (length(intersect(rownames(mat1),rownames(mat2)))==0){stop("The two networks are not interconnected!")}

   # mat_1<-matrix(0,length(matrow),ncol(mat1))
   # rownames(mat_1)<-matrow
   # colnames(mat_1)<-colnames(mat1)
   # mat_1[rownames(mat1),]<-mat1
   #
   # mat_2<-matrix(0,length(matrow),ncol(mat2))
   # rownames(mat_2)<-matrow
   # colnames(mat_2)<-colnames(mat2)
   # mat_2[rownames(mat2),]<-mat2
   if(!weighted) {
   mat1[mat1>0]<-1
   mat2[mat2>0]<-1}
   # mat1<-mat_1
   # mat2<-mat_2

   if(is.null(colnames(mat1)))
      colnames(mat1)<-paste0("a",seq=1:ncol(mat1))
   if(is.null(colnames(mat2)))
      colnames(mat2)<-paste0("c",seq=1:ncol(mat2))
   spe<-unique(c(colnames(mat1),rownames(mat1),rownames(mat2),colnames(mat2)))
   MAT<-matrix(0,length(spe),length(spe))
   dimnames(MAT)<-list(spe,spe)
   MAT[colnames(mat1),rownames(mat1)]<-t(mat1)
   # if(!isDirected1)
   #    MAT[rownames(mat1),colnames(mat1)]<-mat1
   MAT[rownames(mat2),colnames(mat2)]<-mat2
   # if(!isDirected2)
   #    MAT[colnames(mat2),rownames(mat2)]<-t(mat2)
   NET<-graph_from_adjacency_matrix(MAT,weighted=weighted,mode="max")
   V(NET)$name<-spe
   levell<-rep(1,length(spe))
   levell[spe%in%colnames(mat1)]<-0
   # levell[spe%in%rownames(mat1)]<-1
   levell[spe%in%colnames(mat2)]<-2
   V(NET)$level<-levell

   dd<-igraph::layout_with_sugiyama(NET,layers=V(NET)$level)$layout
   dd[order(dd[dd[,2]==3,1]),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==3))
   dd[order(dd[dd[,2]==2,1])+sum(dd[,2]==3),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==2))
   dd[order(dd[dd[,2]==1,1])+sum(dd[,2]==3)+sum(dd[,2]==2),1]<-seq(min(dd[,1]),max(dd[,1]),length.out=sum(dd[,2]==1))
   NET$layout<-dd

   return(NET)
}
