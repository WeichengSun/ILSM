#' Measuring node versatility of multilayer network
#'
#' The versatility of nodes is revealed by several centrality measures that have now been applied to multilayer networks, such as Degree, Pagerank, Hub, Authority, Katz, Eigenvector, and Closeness centrality.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class which contains three groups of species and interactions within layers without interactions between each group of species, or a numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#'
#' @param isDirected1 Logical. Whether the interaction between the two groups of species in \code{mat1} is unidirectional.Default to TRUE, such as Predation and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @param isDirected2 Logical. Whether the interaction between the two groups of species in \code{mat2} is unidirectional.Default to TRUE, such as Predation and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @param type Character. Including "degree", "pagerank", "hub", "authority", "katz", "eigenvector", "closeness", and "all".
#'
#' @details
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
#' \strong{type}
#'
#' \code{type} "degree", "pagerank", "hub", "authority", "katz", "eigenvector", and "closeness" correspond to Degree, PageRank, Hub, Authority, Katz, Eigenvector, and Closeness centrality.
#' \code{type} "all" integrates the above centrality.
#'
#' @return
#'
#' Return a data frame with the first row "node" for each node of network representing each species.
#' \itemize{
#' \item{If \code{type} is either of "degree", "pagerank", "hub", "authority", "katz", "eigenvector", "closeness", the data frame has two columns, and the second column corresponds to either of "Degree", "Pagerank_versatility",
#' "Hub_versatility", "Authority_versatility", "Katz_versatility", "Eigenvector_versatility", "Closeness_versatility".}
#' \item{If \code{type} is "all", the data frame has eight columns, and columns form the second to the eighth correspond to "Degree", "Pagerank_versatility",
#' "Hub_versatility", "Authority_versatility", "Katz_versatility", "Eigenvector_versatility", "Closeness_versatility".}}
#'
#' @references
#' De Domenico, M., Nicosia, V., Arenas, A., & Latora, V. (2015). Structural reducibility of multilayer networks. Nature communications, 6(1), 6864.
#'
#' De Domenico, M., Solé-Ribalta, A., Omodei, E., Gómez, S., & Arenas, A. (2013). Centrality in interconnected multilayer networks. arXiv preprint arXiv:1311.2906.
#'
#' De Domenico, M. (2022). Multilayer Networks: Analysis and Visualization. Introduction to muxViz with R. Cham: Springer.
#'
#' Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The pagerank citation ranking: Bringing order to the web.
#'
#' Magnani, M., Micenkova, B., & Rossi, L. (2013). Combinatorial analysis of multiple networks. arXiv preprint arXiv:1303.4986.
#'
#'
#' @importFrom igraph V
#' @export
#'
#' @examples
#'
#' set.seed(12)
#' d <- build_net(11,22,21,0.2,asmatrices=TRUE)
#' d
#' node_cv(d[[1]])
#'
#' MAT<-d
#' tmat<-t(MAT[[3]])
#' colnames(tmat)<-NULL
#' node_cv(MAT[[3]],MAT[[4]])
#' node_cv(tmat,MAT[[4]])
#' node_cv(MAT[[3]],MAT[[4]],type="pagerank")
#'
#' node_cv(MAT[[3]],MAT[[4]],isDirected2=FALSE)
#'
#'
#'

node_cv <- function(network.or.subnet_mat1,subnet_mat2=NULL,isDirected1=TRUE,isDirected2=TRUE,type=c("degree","pagerank","hub","authority","katz","eigenvector","closeness","all")){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network <- network.or.subnet_mat1
      EDGES <- mEdges(network)
      NODE <- length(network)
      LAYER <- max(EDGES[,4])
      node <- V(network)$name
      MAT<-as.matrix(network[])
      MAT<-MAT+t(MAT)
      MAT[MAT>0]<-1
      MAT<-MAT*upper.tri(MAT)
      degree<-colSums(MAT)+rowSums(MAT)
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
      mat<-edgelist_from_matrices(network.or.subnet_mat1,subnet_mat2,isDirected1=isDirected1,isDirected2=isDirected2)
      EDGES <- mat[[1]]
      node <- mat[[2]]
      NODE <- length(node)
      LAYER <- max(EDGES[,4])
      degree<-mat[[3]]
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")
   SupraAdjacencyMatrix <- BuildSupraAdjacencyMatrixFromExtendedEdgelist_r(mEdges=EDGES, Layers=LAYER, Nodes=NODE, isDirected=F)
   if(missing(type))
      type <- "all"
   if (type == "degree")
      return(data.frame(node=node, degree=degree))
   if (type == "pagerank")
      return(data.frame(node=node, Pagerank_versatility=GetMultiPageRankCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "hub")
      return(data.frame(node=node, Hub_versatility=GetMultiHubCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "authority")
      return(data.frame(node=node, Authority_versatility=GetMultiAuthCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "katz")
      return(data.frame(node=node, Katz_versatility=GetMultiKatzCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "eigenvector")
      return(data.frame(node=node, Eigenvector_versatility=GetMultiEigenvectorCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "closeness")
      return(data.frame(node=node, Closeness_versatility=GetMultiClosenessCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)$closeness))

   if(type == "all"){
      PageRank_versatility <- GetMultiPageRankCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)
      Hub_versatility <- GetMultiHubCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)
      Authority_versatility <- GetMultiAuthCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)
      Katz_versatility <- GetMultiKatzCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)
      Eigenvector_versatility <- GetMultiEigenvectorCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)
      Closeness_versatility <- GetMultiClosenessCentrality_r(SupraAdjacencyMatrix, LAYER, NODE)
      versatility <- data.frame(node=node, Degree=degree, Pagerank_versatility=PageRank_versatility, Hub_versatility=Hub_versatility, Authority_versatility=Authority_versatility,
                                Katz_versatility=Katz_versatility, Eigenvector_versatility=Eigenvector_versatility, Closeness_versatility=Closeness_versatility$closeness)
      rownames(versatility)<-NULL
      return(versatility)
   }
   else
      stop("Error: type is not a valid input!")
}
