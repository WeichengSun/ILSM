#' Measuring Interconnection Centrality of Tripartite Interaction Network
#'
#' The Interconnection Centrality of nodes is revealed by several centrality
#' measures that have now been applied to tripartite networks, such as Degree,
#' Pagerank, Hub, Authority, Eigenvector, and Closeness Betweenness centrality.
#'
#' @param network.or.subnet_mat1 Either a tripartite network of
#' 'igraph' class which contains three groups of species and interactions within
#'  layers without interactions between each group of species, or a numeric
#'  matrix(or data.frame) representing interactions between two groups of
#'  species.
#'  Each row and column of matrix represents single species in the second and
#'  first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are
#'  connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions
#' between two groups of species.
#'  Each row and column of matrix represents single species in the second and
#'  third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are
#'  connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph",
#'  \code{subnet_mat2} defaults to NULL.
#'
#' @param isDirected1 Logical. Whether the interaction between the two groups of
#'  species in \code{mat1} is unidirectional.Default to TRUE, such as Predation
#'  and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @param isDirected2 Logical. Whether the interaction between the two groups of
#'  species in \code{mat2} is unidirectional.Default to TRUE, such as Predation
#'  and Herbivory. Otherwise it is bidirectional, such as Mutualism.
#' @param type Character. Including "degree", "pagerank", "hub", "authority",
#' "eigenvector", "closeness", "betweenness", and "all".
#'
#' @details
#'
#' \strong{network.or.subnet_mat1 and subnet_mat2}
#'
#' There are two types of \code{network.or.subnet_mat1} that can be processed:
#' \itemize{
#' \item Input in a network of type "igraph" alone.
#' \item Must be entered as data frame or matrix with \code{subnet_mat2}.
#' }
#'
#' If the type of inputting is data frame or matrix, please make sure the row of
#'  \code{network.or.subnet_mat1} and \code{subnet_mat2} correspond with the
#'  second group of species that both belong to two subnetworks and interact
#'  with other groups of species.
#' \itemize{
#' \item Try to make the rows of both matrices have the same attributes. Or we
#' default:
#'
#' \item When the two matrices can have different numbers of rows:
#' \itemize{
#' \item 1. If both matrices have row names, then the function counts all row
#' names to produce two new matrices with the same row names.
#' \item 2. If at most one matrix has row names, the function assigns new row
#' names to both matrices on a row-to-row basis (any extra row names are
#' assigned a new value) and then counts all row names to produce two new
#' matrices with the same row names.
#' }
#'
#' \item When the two matrices can have the same numbers of rows:
#' \itemize{
#' \item No matter how the row names of the two matrices are arranged, as long
#' as the row names are exactly the same; But we don't handle matrices with
#' empty row names (the function will give an error).
#' }
#'
#' \item The two matrices can have different numbers of rows, but read our
#' default handling carefully to make sure the calculation is accurate when
#' using this function!!!
#' }
#' About a network of type "igraph", It can be obtained from the connection
#' matrices of subnetworks by the function \code{igraph_from_matrices}.
#'
#' \strong{type}
#'
#' \code{type} "degree", "pagerank", "hub", "authority", "eigenvector",
#' "closeness" and "betweenness" correspond to Degree, PageRank, Hub, Authority,
#'  Eigenvector, Closeness and Betweenness centrality.
#' \code{type} "all" integrates the above centrality.
#'
#' @return
#'
#' Return a data frame with the first row "node" for each node of network
#' representing each species.
#' \itemize{
#' \item If \code{type} is either of "degree", "pagerank", "hub", "authority",
#' "eigenvector", "closeness", "betweenness" the data frame has two columns, and
#'  the second column corresponds to either of "Degree", "Pagerank_Centrality",
#' "Hub_Centrality", "Authority_Centrality", "Eigenvector_Centrality",
#' "Closeness_Centrality", "Betweenness_Centrality".
#' \item If \code{type} is "all", the data frame has eight columns, and columns
#'  form the second to the eighth correspond to "Degree", "Pagerank_Centrality",
#' "Hub_Centrality", "Authority_Centrality", "Eigenvector_Centrality",
#' "Closeness_Centrality", "Betweenness_Centrality".
#' }
#'
#' @references
#' De Domenico, M., Nicosia, V., Arenas, A., & Latora, V. (2015). Structural
#' reducibility of multilayer networks. Nature communications, 6(1), 6864.
#'
#' De Domenico, M., Solé-Ribalta, A., Omodei, E., Gómez, S., & Arenas, A.
#'  (2013). Centrality in interconnected multilayer networks. arXiv preprint
#'  arXiv:1311.2906.
#'
#' De Domenico, M. (2022). Multilayer Networks: Analysis and Visualization.
#' Introduction to muxViz with R. Cham: Springer.
#'
#' Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The pagerank citation
#'  ranking: Bringing order to the web.
#'
#' Magnani, M., Micenkova, B., & Rossi, L. (2013). Combinatorial analysis of
#' multiple networks. arXiv preprint arXiv:1303.4986.
#'
#'
#' @importFrom igraph V
#' @importFrom igraph degree
#' @importFrom igraph page_rank
#' @importFrom igraph hub_score
#' @importFrom igraph authority_score
#' @importFrom igraph centr_eigen
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @export
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' node_icc(Net)
#'
#' MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
#' tmat <- t(MAT[[3]])
#' colnames(tmat) <- NULL
#' node_icc(MAT[[3]],MAT[[4]])
#' node_icc(tmat,MAT[[4]])
#' node_icc(MAT[[3]],MAT[[4]],type="pagerank")
#'
#' node_icc(MAT[[3]],MAT[[4]],isDirected2=FALSE)
#'
#'
#'
#'
node_icc <- function(network.or.subnet_mat1,subnet_mat2=NULL,isDirected1=TRUE,
                     isDirected2=TRUE,type=c("degree","pagerank","hub",
                                             "authority","eigenvector",
                                             "closeness","betweenness","all")){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network <- network.or.subnet_mat1
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T &&
           inherits(subnet_mat2,c("matrix","data.frame"))==T){
      network <- igraph_from_matrices(network.or.subnet_mat1,subnet_mat2,
                                      isDirected1,isDirected2)
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")
   node<-V(network)$name
   Directed<-isDirected1|isDirected2
   if(missing(type))
      type <- "all"
   if (type == "degree")
      return(data.frame(node=node, Degree=degree(network)))
   if (type == "pagerank")
      return(data.frame(node=node, Pagerank_Centrality =
                           page_rank(network, directed = Directed)$vector))
   if (type == "hub")
      return(data.frame(node=node, Hub_Centrality =hub_score(network)$vector))
   if (type == "authority")
      return(data.frame(node=node, Authority_Centrality =
                           authority_score(network)$vector))
   if (type == "eigenvector")
      return(data.frame(node=node, Eigenvector_Centrality =
                           centr_eigen(network)$vector)
             )
   if (type == "closeness")
      return(data.frame(node=node, Closeness_Centrality =
                           closeness(network,mode = "all")))
   if (type == "betweenness")
      return(data.frame(node=node, Betweenness_Centrality =
                           betweenness(network,directed = Directed))
             )

   if(type == "all"){
      Degree <- degree(network)
      PageRank_Centrality  <- page_rank(network,directed = Directed)$vector
      Hub_Centrality  <- hub_score(network)$vector
      Authority_Centrality  <- authority_score(network)$vector
      Eigenvector_Centrality  <- centr_eigen(network)$vector
      Closeness_Centrality  <- closeness(network,mode = "all")
      Betweenness_Centrality  <- betweenness(network,directed = Directed)
      Centrality  <- data.frame(node=node, Degree=Degree, Pagerank_Centrality =
                                   PageRank_Centrality , Hub_Centrality =
                                   Hub_Centrality , Authority_Centrality =
                                   Authority_Centrality ,
                                Eigenvector_Centrality =Eigenvector_Centrality ,
                                Closeness_Centrality =Closeness_Centrality ,
                                Betweenness_Centrality =Betweenness_Centrality )
      rownames(Centrality )<-NULL
      return(Centrality )
   }
   else
      stop("Error: type is not a valid input!")
}
