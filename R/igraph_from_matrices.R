#' Transforming matrices into network
#'
#' Two matrices contain three groups of tropical level species. A multilayer network can be transformed from existing matrices data.
#'
#' @param mat1 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents the species in the first and second layers of the tripartite network respectively.
#'  Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise.
#'
#' @param mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents the species in the second and third layers of the tripartite network respectively.
#'  Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise.
#'
#' @param isDirected1 Logical. Whether the interaction between the two groups of species in \code{mat1} is unidirectional.Default to TRUE, such as Predation and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @param isDirected2 Logical. Whether the interaction between the two groups of species in \code{mat2} is unidirectional.Default to TRUE, such as Predation and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @details
#'
#' \strong{mat1} and \strong{mat2}
#' Make sure that the columns of \code{mat1} correspond strictly to the rows of \code{mat2}:
#' \itemize{
#' \item{The number of columns of \code{mat1} is equal to the number of rows of \code{mat2}.}
#'
#' \item{If both \code{mat1} and \code{mat2} have row and column names, then the column name of \code{mat1} must also correspond to the row name of \code{mat2}.}

#' \item{Or the rows and columns of \code{mat1} and \code{mat2} could be empty.}}
#'
#'
#'
#'
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
#' M <- igraph_from_matrices(MAT[[3]],MAT[[4]])
#' M

igraph_from_matrices<-function(mat1,mat2,isDirected1=TRUE,isDirected2=TRUE){
   if(ncol(mat1)!=nrow(mat2))
      stop("Error: please check whether the column of mat1 is corresponding to the row of mat2!!!")
   if(is.null(colnames(mat1)) || is.null(rownames(mat2))){
      colnames(mat1)<-paste0("mid_pse",seq=1:ncol(mat1))
      rownames(mat2)<-paste0("mid_pse",seq=1:ncol(mat1))
   }
   if(sum(!(colnames(mat1)%in%(rownames(mat2))),na.rm = TRUE)!=0 )
      stop("Error: please check whether the column name of mat1 is corresponding to the row name of mat2!!!")
   if(sum(is.na(colnames(mat1)))!=0 || sum(is.na(rownames(mat2)))!=0)
      stop("Error: There is NA in the column name of mat1 or the row name of mat2!!!")
   mat2<-mat2[colnames(mat1),]
   if(is.null(rownames(mat1)))
      rownames(mat1)<-paste0("up_pse",seq=1:nrow(mat1))
   if(is.null(colnames(mat2)))
      colnames(mat2)<-paste0("down_pse",seq=1:ncol(mat2))
   spe<-unique(c(rownames(mat1),colnames(mat1),rownames(mat2),colnames(mat2)))
   MAT<-matrix(0,length(spe),length(spe))
   dimnames(MAT)<-list(spe,spe)
   MAT[rownames(mat1),colnames(mat1)]<-mat1
   if(!isDirected1)
      MAT[colnames(mat1),rownames(mat1)]<-t(mat1)
   MAT[rownames(mat2),colnames(mat2)]<-mat2
   if(!isDirected2)
      MAT[colnames(mat2),rownames(mat2)]<-t(mat2)
   NET<-graph_from_adjacency_matrix(MAT)
   V(NET)$name<-spe
   levell<-spe
   levell[spe%in%rownames(mat1)]<-0
   levell[spe%in%colnames(mat1)]<-1
   levell[spe%in%colnames(mat2)]<-2
   V(NET)$level<-levell
   return(NET)
}
