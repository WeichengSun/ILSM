#' Intra-guild and inter_guild interaction overlap
#'
#' Calculating species-level intra-guild and inter_guild interaction overlap for a tripartite network with intra-guild interactions.
#'
#' @param mat A square block interaction matrix representing a tripartite network including intra-guild and inter-guild interactions. See details.
#' @param guilds A character vector matching rows of \code{mat} to indicate the guilds using ('a','b' and 'c'). See details.
#' @param method The distance method. Same with \emph{vegan::vegdist}. Default to "horn"
#' @details
#' The input is a block matrix (\eqn{M}) to represent interactions among three groups of species (a-nodes, b-nodes and c-nodes): three intra-guild interaction matrices (\eqn{m_{aa},m_{bb},m_{cc}}),
#' two inter-guild matrices of a and b-nodes (\eqn{m_{ab},m_{ba}} with symmetric links), and two inter-guild matrices of b- and c-nodes(\eqn{m_{bc},m_{cb}} with symmetric links). Connector species belong to b-nodes.
#'
#' \deqn{
#'   \left(
#'     \begin{array}{ccc}
#'       m_{aa} & m_{ab} & 0        \\
#'       m_{ba} & m_{bb} & m_{bc}   \\
#'        0     & m_{cb} & m_{cc}
#'     \end{array}
#'   \right)
#' }
#'
#' \code{guilds} should be a vector of the same length as the row of \code{mat} like c("a","a"..."b","b"..."c","c"..)
#'
#' This function follows the definition by García-Callejas et al (2023). Species-level interaction overlap is derived from the overlap between each pair of species, calculated using the dissimilarity index (\eqn{d_{ij}},default to Morisita–Horn index) as in the R package \emph{vegan}.
#' The net overlap of species (\eqn{o_{i}})is represented by the sum of pairwise overlaps with every other species:
#'   \deqn{
#'          o_i = \sum_{j \in S} (1 - d_{ij})
#'        }

#' @return
#' Return a list including three species-level intra-guild overlap vectors for a-, b- and c-nodes (a_intra_overlap,b_intra_overlap,c_intra_overlap),
#' two vectors of inter-guild interaction overlap for a-nodes (a_inter_b_overlap) and b-nodes (b_inter_a_overlap), and two vectors of inter-guild interaction overlap for b-nodes (b_inter_c_overlap) and c-nodes(c_inter_b_overlap).
#'
#' @references
#' García-Callejas, D., Godoy, O., Buche, L., Hurtado, M., Lanuza, J.B., Allen-Perkins, A. et al. (2023) Non-random interactions within and across guilds shape the potential to coexist in multi-trophic ecological communities. Ecology Letters, 26, 831–842.
#'
#' @examples
#'
#' ## A toy tripartite network with intra-guild negative interactions, inter-guild mutualistic interactions and inter-guild antagonistic interactions.
#' set.seed(12)
#' ##4 a-nodes,5 b-nodes, and 3 c-nodes
#'
#' ##intra-guild interaction matrices
#' mat_aa<-matrix(runif(16,-0.8,-0.2),4,4)
#' mat_bb<-matrix(runif(25,-0.8,-0.2),5,5)
#' mat_cc<-matrix(runif(9,-0.8,-0.2),3,3)
#'
#' ##inter-guild interaction matrices between a- and b-nodes.
#' mat_ab<-mat_ba<-matrix(sample(c(rep(0,8),runif(12,0,0.5))),4,5,byrow=T)# interaction probability = 12/20
#' mat_ba[mat_ba>0]<-runif(12,0,0.5);mat_ba<-t(mat_ba)
#'
#' ##inter-guild interaction matrices between b- and c-nodes.
#' mat_cb<-mat_bc<-matrix(sample(c(rep(0,8),runif(7,0,0.5))),3,5,byrow=T)# interaction probability = 7/15
#' mat_bc[mat_bc>0]<-runif(7,0,0.5);mat_bc<--t(mat_bc)
#' mat<-rbind(cbind(mat_aa,mat_ab,matrix(0,4,3)),cbind(mat_ba,mat_bb,mat_bc),cbind(matrix(0,3,4),mat_cb,mat_cc))
#'
#' ##set the node names
#' rownames(mat)<-c(paste0("a",1:4),paste0("b",1:5),paste0("c",1:3));colnames(mat)<-c(paste0("a",1:4),paste0("b",1:5),paste0("c",1:3))
#' diag(mat)<--1 #assume -1 for diagonal elements
#'
#' ##Visualization of this block matrix.
#' library(plot.matrix)
#' pal <- colorRampPalette(c("darkblue", "lightblue", "white", "pink", "darkred"))(100)
#' par(mar=c(5,5,5,5));plot(mat,col = pal,
#'     breaks = seq(-max(abs(mat)), max(abs(mat)), length.out = 101),  # 101 breaks for 100 colors
#'     main = "Matrix visualization")
#'     clip(x1 = 0.5,# Left boundary
#'     x2 = ncol(mat) + 0.5, # Right boundary
#'     y1 = 0.5,            # Top boundary
#'     y2 = nrow(mat) + 0.5  )
#'     abline(v = c(4.5,9.5), h = c(3.5,8.5), lwd = 3, col = "black")
#'
#' myguilds=c(rep("a",4),rep("b",5),rep("c",3))
#' overlap_guild(mat,guilds=myguilds)

ig_overlap_guild<-function(mat,guilds,method="horn"){

if (!identical(c("a","b","c"),unique(guilds))|length(guilds)!=nrow(mat)){
  Stop("guilds should be a vector including 'a', 'b', 'c'")
}

a_node<-which(guilds=="a")
b_node<-which(guilds=="b")
c_node<-which(guilds=="c")
mat_aa<-mat[a_node,a_node]
mat_bb<-mat[b_node,b_node]
mat_bb<-mat[c_node,c_node]
mat_ab<-mat[a_node,b_node]
mat_ba<-mat[b_node,a_node]
mat_bc<-mat[b_node,c_node]
mat_cb<-mat[c_node,b_node]

a_intra_overlap <- rowSums(1-as.matrix(vegan::vegdist(abs(mat_aa),method = method)),na.rm=T)-1#-1 remove self_overlap
b_intra_overlap <- rowSums(1-as.matrix(vegan::vegdist(abs(mat_bb),method = method)),na.rm=T)-1
c_intra_overlap <- rowSums(1-as.matrix(vegan::vegdist(abs(mat_cc),method = method)),na.rm=T)-1

a_inter_b_overlap<-rowSums(1-as.matrix(vegan::vegdist(abs(mat_ab),method = method)),na.rm=T)
b_inter_a_overlap<-rowSums(1-as.matrix(vegan::vegdist(abs(mat_ba),method = method)),na.rm=T)

b_inter_c_overlap<-rowSums(1-as.matrix(vegan::vegdist(abs(mat_bc),method = method)),na.rm=T)
c_inter_b_overlap<-rowSums(1-as.matrix(vegan::vegdist(abs(mat_cb),method = method)),na.rm=T)

return(list(
       a_intra_overlap=a_intra_overlap,
       b_intra_overlap=b_intra_overlap,
       c_intra_overlap=c_intra_overlap,
       a_inter_b_overlap=a_inter_b_overlap,
       b_inter_a_overlap=b_inter_a_overlap,
       b_inter_c_overlap=b_inter_c_overlap,
       c_inter_b_overlap=c_inter_b_overlap
       ))
}
