#' Transforming matrices into network
#'
#' Two matrices contain three groups of tropical level species. A multilayer network can be transformed from existing matrices data.
#'
#' @param mat1 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param isDirected1 Logical. Whether the interaction between the two groups of species in \code{mat1} is unidirectional.Default to TRUE, such as Predation and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @param isDirected2 Logical. Whether the interaction between the two groups of species in \code{mat2} is unidirectional.Default to TRUE, such as Predation and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @details
#'
#' \strong{mat1} and \strong{mat2}
#'
#' The type of inputting is data frame or matrix, please make sure the row of \code{mat1} and \code{mat2} correspond with the second group of species that both belong to two subnetworks and interact with other groups of species.
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
#'
#' The columns of \code{mat1} and \code{mat2} could be empty. If empty, the function also defaults to the suggested assignment.
#
#' @return
#' Return a network of type "igraph".
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph V
#'
#' @examples
#'
#' set.seed(12)
#' MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
#' MAT[[1]]
#'
#' tmat<-t(MAT[[3]])
#' colnames(tmat)<-NULL
#' igraph_from_matrices(MAT[[3]],MAT[[4]])
#' M <- igraph_from_matrices(tmat,MAT[[4]])
#' M

igraph_from_matrices<-function(mat1, mat2, isDirected1=TRUE, isDirected2=TRUE){
   if(nrow(mat1)!=nrow(mat2))
      message("You will obtain a network that contains non-connecting node!!!")
   else
      message("Because you enter two matrices with the same number of rows, the function will default to them having the same row names. Informed!!")
   if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
      rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
      rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
      matrow<-unique(c(rownames(mat1),rownames(mat2)))
   }
   if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
      matrow<-unique(c(rownames(mat1),rownames(mat2)))
   else
      stop("Make sure matrices either have no row names or have full row names. No NA!!!")
   mat_1<-matrix(0,length(matrow),ncol(mat1))
   rownames(mat_1)<-matrow
   colnames(mat_1)<-colnames(mat1)
   mat_1[rownames(mat1),]<-mat1
   mat_1[mat_1>0]<-1
   mat_2<-matrix(0,length(matrow),ncol(mat2))
   rownames(mat_2)<-matrow
   colnames(mat_2)<-colnames(mat2)
   mat_2[rownames(mat2),]<-mat2
   mat_2[mat_2>0]<-1
   mat1<-mat_1
   mat2<-mat_2

   if(is.null(colnames(mat1)))
      colnames(mat1)<-paste0("up_spe",seq=1:ncol(mat1))
   if(is.null(colnames(mat2)))
      colnames(mat2)<-paste0("down_spe",seq=1:ncol(mat2))
   spe<-unique(c(colnames(mat1),rownames(mat1),rownames(mat2),colnames(mat2)))
   MAT<-matrix(0,length(spe),length(spe))
   dimnames(MAT)<-list(spe,spe)
   MAT[colnames(mat1),rownames(mat1)]<-t(mat1)
   if(!isDirected1)
      MAT[rownames(mat1),colnames(mat1)]<-mat1
   MAT[rownames(mat2),colnames(mat2)]<-mat2
   if(!isDirected2)
      MAT[colnames(mat2),rownames(mat2)]<-t(mat2)
   NET<-graph_from_adjacency_matrix(MAT)
   V(NET)$name<-spe
   levell<-spe
   levell[spe%in%colnames(mat1)]<-0
   levell[spe%in%rownames(mat1)]<-1
   levell[spe%in%colnames(mat2)]<-2
   V(NET)$level<-levell
   return(NET)
}
