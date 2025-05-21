#' Interconnection centrality for connector nodes in a tripartite network
#'
#' Calculating interconnection centrality (degree, betweenness and closeness) specified for connector nodes in a tripartite network.
#'
#' @param network.or.subnet_mat1 Igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork.
#' @param subnet_mat2 The matrix representing one subnetwork.
#' @param weighted Logical. If TRUE, link strengths of connector nodes are used. Default to FALSE.

#' #' @details
#'This function calculates interconnection degree, betweenness and closeness centrality for connector nodes.
#'For binary networks, interconnection degree centrality of each connector species is defined as the product of its degree values from two subnetworks,
#'interconnection betweenness centrality is defined by the number of shortest paths from a-nodes to c-nodes going through connector species,
#'and interconnection closeness centrality is defined by the inverse of the sum of distances from connector species to both a-nodes and c-nodes.
#'For weighted networks, interaction strengths are used in the calculation of weighted degree, shorest path, and distance.
#' \strong{network.or.subnet_mat1 and subnet_mat2}
#'
#' \strong{weighted}
#'
#' If the \code{weighted} = FALSE, the input network can be an "igraph" object or two matrices.If a weighted network is provided, it will be transformed to a binary network.
#' If the \code{weighted} = TRUE, the input network can only be two matrices.
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
#'
#' Return a data frame with interconnection degree, betweenness and closeness.
#' @references
#' De Domenico, M., Nicosia, V., Arenas, A., & Latora, V. (2015). Structural
#' reducibility of multilayer networks. Nature communications, 6(1), 6864.
#'
#' De Domenico, M., SolC)-Ribalta, A., Omodei, E., GC3mez, S., & Arenas, A.
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
#' set.seed(13)
#' E(Net)$weight<-runif(length(E(Net)),0.1,1)#random weights assigned
#' node_icc(Net,weighted=T)
#'
#'##input as binary matrices,with row names.
#' set.seed(12)
#' md1 <- matrix(sample(c(0,1),5*8,replace=TRUE),5,8,dimnames = list(paste0("b",1:5),paste0("c",1:8)))
#' md2 <- matrix(sample(c(0,1),20*30,replace=TRUE),20,30,dimnames = list(paste0("b",1:20),paste0("a",1:30)))
#' node_icc(md1,md2)

#'
#'##input as weighted matrices,with row numbers as row names.
#' set.seed(17)
#' mdw1 <- matrix(sample(c(rep(0,20),runif(20,0,1))),5,8)
#' mdw2 <- matrix(sample(c(rep(0,500),runif(100,0,1))),20,30)
#' node_icc(mdw1,mdw2)
#' node_icc(mdw1,mdw2,weighted=T)
#'
node_icc <- function(network.or.subnet_mat1,subnet_mat2=NULL,weighted=F){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network <- network.or.subnet_mat1
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T &&
           inherits(subnet_mat2,c("matrix","data.frame"))==T){
      network <- igraph_from_matrices(network.or.subnet_mat1,subnet_mat2,weighted=weighted)
   }
   else
      stop("please check the type of 'network.or.subnet_mat1'")

   network0 <- adjust_net(network,weighted=T)
   connector_node=V(network0)$name[V(network0)$level==1]
   mat <- as.matrix(network[])

   mat1 <- mat[V(network)$name[V(network)$level%in% (0:1)],V(network)$name[V(network)$level%in% (0:1)]]
   net1 <- graph_from_adjacency_matrix(mat1,weighted = T,mode = "max")


   mat2 <- mat[V(network)$name[V(network)$level%in%(1:2)],V(network)$name[V(network)$level %in% (1:2)]]
   net2 <- graph_from_adjacency_matrix(mat2,weighted = T,mode = "max")

   if(!weighted){
      E(net1)$weight=1
      E(net2)$weight=1
   }
   net1_degree <- strength(net1)[connector_node]
   net2_degree <- strength(net2)[connector_node]
   net_degree <- net1_degree * net2_degree

   if(!weighted){
      E(network)$weight=1
      }
   #net_closeness
   dist_res<-distances(network,v=connector_node,to=V(network)$name[V(network)$level%in% c(0,2)],mode="all")
   dist_res<-as.data.frame(dist_res)
   net_closeness<-1/apply(dist_res,1,function(x){sum(x[!is.na(x)&!is.infinite(x)])})##有INF列

   #net_betweenness

   #tmp<-lapply(V(network)$name[V(network)$level%in% c(0)],function(i){all_shortest_paths(network,from=i,to=V(network)$name[V(network)$level%in% c(2)])$vpaths})
   #tmp<-unlist(tmp,recursive = F)
   connector_ID=as.vector(V(network)[V(network)$name%in%connector_node])
   b<-as.vector(V(network)[V(network)$level%in% c(1)])
   a<-as.vector(V(network)[V(network)$level%in% c(0)])
   c<-as.vector(V(network)[V(network)$level%in% c(2)])
   tmp<-apply(expand.grid(a,c),1,function(p){do.call(rbind,all_shortest_paths(network,from=p[1],to=p[2])$vpaths)})
   tmp<-tmp[!sapply(tmp,is.null)]

   res<-sapply(tmp,function(x){
      x<-x[,2:(ncol(x)-1),drop=F]
      res0<-apply(x,1,function(y) connector_ID%in%y)
      rowSums(res0)/nrow(x)

   })
   #tmp<-data.frame(do.call(rbind,tmp))
   #names(tmp)<-c("a","b","c")

   # res<-apply(expand.grid(a,c),1,function(p){
   #               is_ac<-sapply(tmp,function(x){x[1]==p[1]&x[length(x)]==p[2]})
   #               g_ac<-sum(is_ac)#nrow(tmp[tmp$a==p[1]&tmp$c==p[2],])
   #               sapply(b,function(x){sapply(tmp[is_ac], function(y) {y[-c(1,length(y))]%in%x})/g_ac})
   #            })

   net_betweenness<-rowSums(res,na.rm=T)

   Centrality  <- data.frame(node=connector_node, interconnection_degree=net_degree, interconnection_betweenness=net_betweenness,interconnection_closeness=net_closeness)
   return(Centrality)
}
