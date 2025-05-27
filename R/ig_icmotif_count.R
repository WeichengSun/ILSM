#' Count interconnection motifs for tripartite networks with intra-guild interactions
#'
#' Counting the frequencies of interconnection motifs for a tripartite interaction network with intra-guild interactions.
#'
#' @param network.or.subnet_mat1 An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.
#' @param subnet_mat2 A matrix representing one subnetwork.
#' @param weighted Logical. Default to FALSE. If TRUE, the arithmetic mean of the subgraph weights is provided for each motif. See details
#' @import igraph
#'
#' @export
#'
#' @details
#' This function is designed for tripartite networks with intra-guild interactions. The input network should be nput as a block matrix (\eqn{M}) to represent three groups of nodes (a-nodes, b-nodes and c-nodes): three intra-guild interaction matrices (\eqn{m_{aa},m_{bb},m_{cc}}),
#' two inter-guild matrices of a and b-nodes (\eqn{m_{ab},m_{ba}}), and two inter-guild matrices of b- and c-nodes(\eqn{m_{bc},m_{cb}}).
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
#' \strong{Interconnection motifs in tripartite networks with intra-guild interactions}
#' <br>An interconnection motif is defined to comprise three sets of connected nodes: the connector nodes (belonging to b-nodes), the nodes in one subnetwork (belonging to a-nodes in the P subnetwork), and the nodes in the other subnetwork (belonging to c-nodes in the Q subnetwork). Each guild has two nodes at most, resulting in a total of 107 distinct motif forms.
#' The algorithm for counting interconnection motifs is designed by extending the fast approach from Simmons et al.(2019). For interconnection motifs in tripartite networks without intra-guild interactions, please see **icmotif_count** and **icmotif_role**.
#'
#' \strong{Weighted networks}
#' <br>For weighted tripartite networks, the mean weight of a given motif is provided by averaging the weights of all motif occurrences isomorphic to the motif. The weight of a motif occurrence is the arithmetic mean of the weights of its links, following Mora et al. (2018) and Simmons et al. (2019).
#'
#'
#' @return
#'  Return a data.fame of the frequencies (and mean weight) of 107 interconnection motifs.
#'
#' @srrstats {G1.1} The algorithm is the first implementation of a novel algorithm.
#' @srrstats {G1.3,G1.4} This standard belongs here.
#' @srrstats {G2.1} It can input single- or multi- data.
#' @srrstats {G2.1a} Because the input data flexibility of \code{network.or.subnet_mat1} and \code{subnet_mat2}, it provide explicit secondary documentation.
#' @srrstats {G2.4,G2.4a} This standard belongs here.
#' @srrstats {G2.7} It accepts two types of \code{network.or.subnet_mat1}.
#' @srrstats {G2.13} It checks the possibility of ('NA') data.
#' @srrstats {G2.14,G2.14a} For the missing ('NA') data, it provide the error on missing data.
#' @srrstats {G2.15} This function never assume non-missingness, and never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters.
#' @srrstats {G3.0} All numeric equality comparisons ensure that they are made between integers
#' @srrstats {G5.0} The software contains standard data sets with known properties, such as \code{data: PPH_Coltparkmeadow}.
#' @srrstats {G5.2,G5.2a} Here the message produced within R code by 'stop()'.
#'
#' @references
#' Pilosof, S., Porter, M. A., Pascual, M., & KC)fi, S. (2017). The multilayer nature of ecological networks. Nature Ecology & Evolution, 1(4), 0101.
#'
#' Mora, B.B., Cirtwill, A.R. and Stouffer, D.B. (2018). pymfinder: a tool for the motif analysis of binary and quantitative complex networks. bioRxiv, 364703.
#'
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#' icmotif_count(Net)
#'
#' ## empirical network
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' icmotif_count(Net)
#' set.seed(13)
#' E(Net)$weight<-runif(length(E(Net)),0.1,1)#random weights assigned
#' icmotif_count(Net, weighted=T)
#'
#'
#' ##input as binary matrices, with row names.
#' set.seed(12)
#' md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
#' md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
#' icmotif_count(md1,md2)
#'
#'
#'##input as weighted matrices,with row numbers as row names.
#' set.seed(12)
#' mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
#' mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' icmotif_count(mdw1,mdw2,weighted=T)
#'
ig_icmotif_count<-function(subnet_mat1, subnet_mat2,
            AA= NULL, BB= NULL, CC= NULL,
            weighted =FALSE) {
      if (inherits(subnet_mat1, c("matrix", "data.frame")) == T &&
          inherits(subnet_mat2, c("matrix", "data.frame")) == T ){
         if(is.null(AA)){
            AA <- matrix(0, nrow = ncol(subnet_mat1), ncol = ncol(subnet_mat1))
            dimnames(AA) <- list(colnames(subnet_mat1),colnames(subnet_mat1))
         }
         if(is.null(BB)){
            BB <- matrix(0, nrow = nrow(subnet_mat1), ncol = nrow(subnet_mat1))
            dimnames(BB) <- list(rownames(subnet_mat1),rownames(subnet_mat1))
         }
         if(is.null(CC)){
            CC <- matrix(0, nrow = ncol(subnet_mat2), ncol = ncol(subnet_mat2))
            dimnames(CC) <- list(colnames(subnet_mat2),colnames(subnet_mat2))
         }

         if (inherits(AA, c("matrix", "data.frame")) == T &&
             inherits(BB, c("matrix", "data.frame")) == T &&
             inherits(CC, c("matrix", "data.frame")) == T) {
            if (is.null(rownames(subnet_mat1)) |
                is.null(rownames(subnet_mat2)) | is.null(rownames(BB)))
               stop("Make sure these matrices either have no row names or have full row names. [No NA!!!]")
            matrow <- unique(c(rownames(subnet_mat1), rownames(subnet_mat2)))

            matrow_l <- length(matrow)
            Brow <- rownames(BB)[rownames(BB) %in% matrow]
            BB_mat <- matrix(0, matrow_l, matrow_l)
            dimnames(BB_mat) <- list(matrow, matrow)
            BB_mat[Brow, Brow] <- BB[Brow, Brow]


            if (is.null(colnames(subnet_mat1)) | is.null(rownames(AA)))
               stop("Make sure these matrices either have no row and column names or have full row and column names. [No NA!!!]")
            A_l <- ncol(subnet_mat1)
            Arow <- rownames(AA)[rownames(AA) %in% colnames(subnet_mat1)]
            AA_mat <- matrix(0, A_l, A_l)
            dimnames(AA_mat) <-
               list(colnames(subnet_mat1), colnames(subnet_mat1))
            AA_mat[Arow, Arow] <- AA[Arow, Arow]

            if (is.null(colnames(subnet_mat2)) | is.null(rownames(CC)))
               stop("Make sure these matrices either have no row and column names or have full row and columnnames. [No NA!!!]")
            C_l <- ncol(subnet_mat2)
            Crow <- rownames(CC)[rownames(CC) %in% colnames(subnet_mat2)]
            CC_mat <- matrix(0, C_l, C_l)
            dimnames(CC_mat) <-
               list(colnames(subnet_mat2), colnames(subnet_mat2))
            CC_mat[Crow, Crow] <- CC[Crow, Crow]
         }
      }

      mat_1 <- matrix(0, matrow_l, A_l)
      rownames(mat_1) <- matrow
      mat_1[rownames(subnet_mat1), ] <- subnet_mat1
      mat_2 <- matrix(0, matrow_l, C_l)
      rownames(mat_2) <- matrow
      mat_2[rownames(subnet_mat2), ] <- subnet_mat2
      dimnames(mat_1) <- NULL
      dimnames(mat_2) <- NULL
      P <- t(mat_1)
      Q <- mat_2
      AA <- AA_mat
      BB <- BB_mat
      CC <- CC_mat

      PW <- P; QW <- Q; AAW <- AA; BBW <- BB; CCW <- CC
      AA <- (AA!=0)*1; BB <- (BB!=0)*1; CC <- (CC!=0)*1
      P <- (P!=0)*1; Q <- (Q!=0)*1

      AA_no <- AA == 0
      diag(AA) <- 0
      diag(AA_no) <- 0
      AA_two <- (AA!=0)*1

      BB_no <- BB == 0
      diag(BB) <- 0
      diag(BB_no) <- 0
      BB_two <- (BB!=0)*1


      CC_no <- CC == 0
      diag(CC) <- 0
      diag(CC_no) <- 0
      CC_two <- (CC!=0)*1


      Ob <- colSums(P); Rb <- rowSums(Q)
      Ubb <- t(P) %*% P
      Vbb <- Q %*% t(Q)
      UVbb <- Ubb * Vbb
      B2C_NO <- rowSums((Q%*%CC_no) *Q) /2
      B2C_YES <- rowSums((Q%*%CC_two) *Q) /2
      B2A_NO <- rowSums((t(P)%*%AA_no) *t(P)) /2
      B2A_YES <- rowSums((t(P)%*%AA_two) *t(P)) /2


      M111 <- sum(Ob * Rb)
      M112_1 <- sum(Ob * B2C_NO)
      M112_2 <- sum(Ob * B2C_YES)
      M211_1 <- sum(B2A_NO * Rb)
      M211_2 <- sum(B2A_YES * Rb)
      M212_1 <- sum(B2A_NO * B2C_NO)
      M212_2 <- sum(B2A_YES * B2C_NO)
      M212_3 <- sum(B2A_NO * B2C_YES)
      M212_4 <- sum(B2A_YES * B2C_YES)

      mat_tri <- function(mat) {sum(mat[upper.tri(mat)])}

      M121_1 <- mat_tri(UVbb * BB_no)
      M121_2 <- mat_tri(UVbb * BB_two)



      P_all_yes <- t(apply(P, 2, function(x) { A <- x %*% t(x); A * AA}))
      P_one_yes <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); A *AA}))
      P_two_yes <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); t(A) *AA}))

      P_all_no <- t(apply(P, 2, function(x) { A <- x %*% t(x); A * AA_no}))
      P_one_no <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); A * AA_no}))
      P_two_no <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); t(A) * AA_no}))

      Q_all_yes <- t(apply(Q, 1, function(x) { A <- x %*% t(x); A * CC}))
      Q_one_yes <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); A * CC}))
      Q_two_yes <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); t(A) * CC}))

      Q_all_no <- t(apply(Q, 1, function(x) { A <- x %*% t(x); A * CC_no}))
      Q_one_no <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); A * CC_no}))
      Q_two_no <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); t(A) * CC_no}))

      AA_no_P_2 <- P_two_no %*% t(P_one_no)
      AA_yes_P_2 <- P_two_yes %*% t(P_one_yes)
      AA_no_P_3 <- P_all_no %*% t(P_one_no)
      AA_yes_P_3 <- P_all_yes %*% t(P_one_yes)
      AA_no_P_3t <- t(AA_no_P_3)
      AA_yes_P_3t <- t(AA_yes_P_3)
      AA_no_P_4 <- (P_all_no %*% t(P_all_no)) / 2
      AA_yes_P_4 <- (P_all_yes %*% t(P_all_yes)) / 2


      CC_no_Q_2 <- Q_two_no %*% t(Q_one_no)
      CC_yes_Q_2 <- Q_two_yes %*% t(Q_one_yes)
      CC_no_Q_3 <- Q_all_no %*% t(Q_one_no)
      CC_yes_Q_3 <- Q_all_yes %*% t(Q_one_yes)
      CC_no_Q_4 <- (Q_all_no %*% t(Q_all_no)) / 2
      CC_yes_Q_4 <- (Q_all_yes %*% t(Q_all_yes)) / 2



      M122_1_1 <- mat_tri(Ubb * BB_no * CC_no_Q_2)
      M122_1_2 <- mat_tri(Ubb * BB * CC_no_Q_2)
      M122_1_3 <- mat_tri(Ubb * BB_no * CC_yes_Q_2)
      M122_1_4 <- mat_tri(Ubb * BB * CC_yes_Q_2)
      M122_2_1 <- sum(Ubb * BB_no * CC_no_Q_3)
      M122_2_2 <- sum(Ubb * BB * CC_no_Q_3)
      M122_2_3 <- sum(Ubb * BB_no * CC_yes_Q_3)
      M122_2_4 <- sum(Ubb * BB * CC_yes_Q_3)
      M122_3_1 <- mat_tri(Ubb * BB_no * CC_no_Q_4)
      M122_3_2 <- mat_tri(Ubb * BB * CC_no_Q_4)
      M122_3_3 <- mat_tri(Ubb * BB_no * CC_yes_Q_4)
      M122_3_4 <- mat_tri(Ubb * BB * CC_yes_Q_4)


      M221_1_1 <- mat_tri(AA_no_P_2 * BB_no * Vbb)
      M221_1_2 <- mat_tri(AA_yes_P_2 * BB_no * Vbb)
      M221_1_3 <- mat_tri(AA_no_P_2 * BB * Vbb)
      M221_1_4 <- mat_tri(AA_yes_P_2 * BB * Vbb)
      M221_2_1 <- sum(AA_no_P_3 * BB_no * Vbb)
      M221_2_2 <- sum(AA_yes_P_3 * BB_no * Vbb)
      M221_2_3 <- sum(AA_no_P_3 * BB * Vbb)
      M221_2_4 <- sum(AA_yes_P_3 * BB * Vbb)
      M221_3_1 <- mat_tri(AA_no_P_4 * BB_no * Vbb)
      M221_3_2 <- mat_tri(AA_yes_P_4 * BB_no * Vbb)
      M221_3_3 <- mat_tri(AA_no_P_4 * BB * Vbb)
      M221_3_4 <- mat_tri(AA_yes_P_4 * BB * Vbb)

      M222_1_1 <- sum(AA_no_P_3 * BB_no * CC_no_Q_2)
      M222_1_2 <- sum(AA_yes_P_3 * BB_no * CC_no_Q_2)
      M222_1_3 <- sum(AA_no_P_3 * BB * CC_no_Q_2)
      M222_1_4 <- sum(AA_no_P_3 * BB_no * CC_yes_Q_2)
      M222_1_5 <- sum(AA_yes_P_3 * BB * CC_no_Q_2)
      M222_1_6 <- sum(AA_yes_P_3 * BB_no * CC_yes_Q_2)
      M222_1_7 <- sum(AA_no_P_3 * BB * CC_yes_Q_2)
      M222_1_8 <- sum(AA_yes_P_3 * BB * CC_yes_Q_2)

      M222_2_1 <- mat_tri(AA_no_P_4 * BB_no * CC_no_Q_2)
      M222_2_2 <- mat_tri(AA_yes_P_4 * BB_no * CC_no_Q_2)
      M222_2_3 <- mat_tri(AA_no_P_4 * BB * CC_no_Q_2)
      M222_2_4 <- mat_tri(AA_no_P_4 * BB_no * CC_yes_Q_2)
      M222_2_5 <- mat_tri(AA_yes_P_4 * BB * CC_no_Q_2)
      M222_2_6 <- mat_tri(AA_yes_P_4 * BB_no * CC_yes_Q_2)
      M222_2_7 <- mat_tri(AA_no_P_4 * BB * CC_yes_Q_2)
      M222_2_8 <- mat_tri(AA_yes_P_4 * BB * CC_yes_Q_2)

      M222_3_1 <- sum(AA_no_P_2 * BB_no * CC_no_Q_3)
      M222_3_2 <- sum(AA_yes_P_2 * BB_no * CC_no_Q_3)
      M222_3_3 <- sum(AA_no_P_2 * BB * CC_no_Q_3)
      M222_3_4 <- sum(AA_no_P_2 * BB_no * CC_yes_Q_3)
      M222_3_5 <- sum(AA_yes_P_2 * BB * CC_no_Q_3)
      M222_3_6 <- sum(AA_yes_P_2 * BB_no * CC_yes_Q_3)
      M222_3_7 <- sum(AA_no_P_2 * BB * CC_yes_Q_3)
      M222_3_8 <- sum(AA_yes_P_2 * BB * CC_yes_Q_3)



      M222_4_1 <- sum(AA_no_P_3 * BB_no * CC_no_Q_3)
      M222_4_2 <- sum(AA_yes_P_3 * BB_no * CC_no_Q_3)
      M222_4_3 <- sum(AA_no_P_3 * BB * CC_no_Q_3)
      M222_4_4 <- sum(AA_no_P_3 * BB_no * CC_yes_Q_3)
      M222_4_5 <- sum(AA_yes_P_3 * BB * CC_no_Q_3)
      M222_4_6 <- sum(AA_yes_P_3 * BB_no * CC_yes_Q_3)
      M222_4_7 <- sum(AA_no_P_3 * BB * CC_yes_Q_3)
      M222_4_8 <- sum(AA_yes_P_3 * BB * CC_yes_Q_3)



      M222_5_1 <- sum(AA_no_P_4 * BB_no * CC_no_Q_3)
      M222_5_2 <- sum(AA_yes_P_4 * BB_no * CC_no_Q_3)
      M222_5_3 <- sum(AA_no_P_4 * BB * CC_no_Q_3)
      M222_5_4 <- sum(AA_no_P_4 * BB_no * CC_yes_Q_3)
      M222_5_5 <- sum(AA_yes_P_4 * BB * CC_no_Q_3)
      M222_5_6 <- sum(AA_yes_P_4 * BB_no * CC_yes_Q_3)
      M222_5_7 <- sum(AA_no_P_4 * BB * CC_yes_Q_3)
      M222_5_8 <- sum(AA_yes_P_4 * BB * CC_yes_Q_3)



      M222_6_1 <- mat_tri(AA_no_P_2 * BB_no * CC_no_Q_4)
      M222_6_2 <- mat_tri(AA_yes_P_2 * BB_no * CC_no_Q_4)
      M222_6_3 <- mat_tri(AA_no_P_2 * BB * CC_no_Q_4)
      M222_6_4 <- mat_tri(AA_no_P_2 * BB_no * CC_yes_Q_4)
      M222_6_5 <- mat_tri(AA_yes_P_2 * BB * CC_no_Q_4)
      M222_6_6 <- mat_tri(AA_yes_P_2 * BB_no * CC_yes_Q_4)
      M222_6_7 <- mat_tri(AA_no_P_2 * BB * CC_yes_Q_4)
      M222_6_8 <- mat_tri(AA_yes_P_2 * BB * CC_yes_Q_4)


      M222_7_1 <- sum(AA_no_P_3 * BB_no * CC_no_Q_4)
      M222_7_2 <- sum(AA_yes_P_3 * BB_no * CC_no_Q_4)
      M222_7_3 <- sum(AA_no_P_3 * BB * CC_no_Q_4)
      M222_7_4 <- sum(AA_no_P_3 * BB_no * CC_yes_Q_4)
      M222_7_5 <- sum(AA_yes_P_3 * BB * CC_no_Q_4)
      M222_7_6 <- sum(AA_yes_P_3 * BB_no * CC_yes_Q_4)
      M222_7_7 <- sum(AA_no_P_3 * BB * CC_yes_Q_4)
      M222_7_8 <- sum(AA_yes_P_3 * BB * CC_yes_Q_4)

      M222_8_1 <- mat_tri(AA_no_P_4 * BB_no * CC_no_Q_4)
      M222_8_2 <- mat_tri(AA_yes_P_4 * BB_no * CC_no_Q_4)
      M222_8_3 <- mat_tri(AA_no_P_4 * BB * CC_no_Q_4)
      M222_8_4 <- mat_tri(AA_no_P_4 * BB_no * CC_yes_Q_4)
      M222_8_5 <- mat_tri(AA_yes_P_4 * BB * CC_no_Q_4)
      M222_8_6 <- mat_tri(AA_yes_P_4 * BB_no * CC_yes_Q_4)
      M222_8_7 <- mat_tri(AA_no_P_4 * BB * CC_yes_Q_4)
      M222_8_8 <- mat_tri(AA_yes_P_4 * BB * CC_yes_Q_4)



      M222_9_1 <- sum(AA_no_P_3t * BB_no * CC_no_Q_3)
      M222_9_2 <- sum(AA_yes_P_3t * BB_no * CC_no_Q_3)
      M222_9_3 <- sum(AA_no_P_3t * BB * CC_no_Q_3)
      M222_9_4 <- sum(AA_no_P_3t * BB_no * CC_yes_Q_3)
      M222_9_5 <- sum(AA_yes_P_3t * BB * CC_no_Q_3)
      M222_9_6 <- sum(AA_yes_P_3t * BB_no * CC_yes_Q_3)
      M222_9_7 <- sum(AA_no_P_3t * BB * CC_yes_Q_3)
      M222_9_8 <- sum(AA_yes_P_3t * BB * CC_yes_Q_3)

      intra_guild_motif <- c(M111,M112_1,M112_2,M211_1,M211_2,M212_1,M212_2,M212_3,M212_4,M121_1,M121_2,M122_1_1,M122_1_2,M122_1_3,M122_1_4,
               M122_2_1,M122_2_2,M122_2_3,M122_2_4,M122_3_1,M122_3_2,M122_3_3,M122_3_4,M221_1_1,M221_1_2,M221_1_3,M221_1_4,
               M221_2_1,M221_2_2,M221_2_3,M221_2_4,M221_3_1,M221_3_2,M221_3_3,M221_3_4,M222_1_1,M222_1_2,M222_1_3,M222_1_4,
               M222_1_5,M222_1_6,M222_1_7,M222_1_8,M222_2_1,M222_2_2,M222_2_3,M222_2_4,M222_2_5,M222_2_6,M222_2_7,M222_2_8,
               M222_3_1,M222_3_2,M222_3_3,M222_3_4,M222_3_5,M222_3_6,M222_3_7,M222_3_8,
               M222_4_1,M222_4_2,M222_4_3,M222_4_4,M222_4_5,M222_4_6,M222_4_7,M222_4_8,
               M222_5_1,M222_5_2,M222_5_3,M222_5_4,M222_5_5,M222_5_6,M222_5_7,M222_5_8,
               M222_6_1,M222_6_2,M222_6_3,M222_6_4,M222_6_5,M222_6_6,M222_6_7,M222_6_8,
               M222_7_1,M222_7_2,M222_7_3,M222_7_4,M222_7_5,M222_7_6,M222_7_7,M222_7_8,
               M222_8_1,M222_8_2,M222_8_3,M222_8_4,M222_8_5,M222_8_6,M222_8_7,M222_8_8,
               M222_9_1,M222_9_2,M222_9_3,M222_9_4,M222_9_5,M222_9_6,M222_9_7,M222_9_8)

      motif_names<-c("M111","M112-1","M112-2","M211-1","M211-2","M212-1","M212-2","M212-3","M212-4","M121-1","M121-2","M122-1-1","M122-1-2","M122-1-3","M122-1-4",
                     "M122-2-1","M122-2-2","M122-2-3","M122-2-4","M122-3-1","M122-3-2","M122-3-3","M122-3-4","M221-1-1","M221-1-2","M221-1-3","M221-1-4",
                     "M221-2-1","M221-2-2","M221-2-3","M221-2-4","M221-3-1","M221-3-2","M221-3-3","M221-3-4","M222-1-1","M222-1-2",'M222-1-3',"M222-1-4",
                     "M222-1-5","M222-1-6","M222-1-7","M222-1-8","M222-2-1","M222-2-2","M222-2-3","M222-2-4","M222-2-5","M222-2-6","M222-2-7","M222-2-8",
                     "M222-3-1","M222-3-2","M222-3-3","M222-3-4","M222-3-5","M222-3-6","M222-3-7","M222-3-8",
                     "M222-4-1","M222-4-2","M222-4-3","M222-4-4","M222-4-5","M222-4-6","M222-4-7","M222-4-8",
                     "M222-5-1","M222-5-2","M222-5-3","M222-5-4","M222-5-5","M222-5-6","M222-5-7","M222-5-8",
                     "M222-6-1","M222-6-2","M222-6-3","M222-6-4","M222-6-5","M222-6-6","M222-6-7","M222-6-8",
                     "M222-7-1","M222-7-2","M222-7-3","M222-7-4","M222-7-5","M222-7-6","M222-7-7","M222-7-8",
                     "M222-8-1","M222-8-2","M222-8-3","M222-8-4","M222-8-5","M222-8-6","M222-8-7","M222-8-8",
                     "M222-9-1","M222-9-2","M222-9-3","M222-9-4","M222-9-5","M222-9-6","M222-9-7","M222-9-8")
      if(weighted){

         All_add <- function(a,b){ return((a!=0)*(b!=0)*(a+b)) }
         One_add <- function(a,b){ return((a!=0)*(b!=0)*(a)) }

         ObW <- colSums(PW)
         RbW <- rowSums(QW)

         B2C_NOW <- apply(QW,1,function(x){ A<- apply(CC_no,2,function(y){ One_add(x,y) });  return(sum(All_add(x,t(A)))/2)})
         B2C_YESW <- apply(QW,1,function(x){ A<- apply(CCW,2,function(y){ All_add(x,y) });  return(sum(All_add(x,t(A)))/2)})

         B2A_NOW <- apply(PW,2,function(x){ A<- apply(AA_no,2,function(y){ One_add(x,y) });  return(sum(All_add(x,t(A)))/2)})
         B2A_YESW <- apply(PW,2,function(x){ A<- apply(AAW,2,function(y){ All_add(x,y) });  return(sum(All_add(x,t(A)))/2)})


         M111 <- sum(ObW * Rb + RbW * Ob)

         M112_1 <- sum( ObW * B2C_NO + B2C_NOW * Ob  )
         M112_2 <- sum( ObW * B2C_YES + B2C_YESW * Ob  )

         M211_1 <- sum(B2A_NOW * Rb + RbW * B2A_NO )
         M211_2 <- sum(B2A_YESW * Rb +  RbW * B2A_YES )

         M212_1 <- sum(B2A_NOW * B2C_NO + B2C_NOW *B2A_NO)
         M212_2 <- sum(B2A_YESW * B2C_NO + B2C_NOW * B2A_YES)
         M212_3 <- sum(B2A_NOW * B2C_YES + B2C_YESW * B2A_NO)
         M212_4 <- sum(B2A_YESW * B2C_YES + B2C_YESW * B2A_YES)


         A_2BW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum(All_add(x,y)) })}))
         C_2BW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum(All_add(x,y)) })}))


         M121_1 <- mat_tri( (A_2BW*Vbb + C_2BW*Ubb) * BB_no)
         M121_2 <- mat_tri( (A_2BW*Vbb + C_2BW*Ubb)* BB_two + UVbb * BBW )



         PW_all_yes <- t(apply(PW, 2, function(x){A <- (x!=0)%*%t(x); All_add(A,x)}))  # All_add(t(All_add(t(All_add(A,x)),x)),PP)
         PW_one_yes <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); A }))
         PW_two_yes <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); t(A) }))

         PW_all_no <- t(apply(PW, 2, function(x){A <- (x!=0)%*%t(x);  All_add(A,x) * AA_no}))  # All_add(t(All_add(t(All_add(A,x)),x)),PP)
         PW_one_no <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); A * AA_no}))
         PW_two_no <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); t(A) * AA_no}))

         QW_all_yes <- t(apply(QW, 1, function(x){A <- (x!=0)%*%t(x); All_add(A,x)}))
         QW_one_yes <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); A}))
         QW_two_yes <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); t(A) }))

         QW_all_no <- t(apply(QW, 1, function(x){A <- (x!=0)%*%t(x); All_add(A,x)* CC_no}))
         QW_one_no <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); A * CC_no}))
         QW_two_no <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); t(A) * CC_no}))




         AA_no_PW_2 <- t(apply(PW_two_no,1,function(x){apply(PW_one_no,1,function(y){ sum(All_add(x,y))})}))
         AA_yes_PW_2 <- t(apply(PW_two_yes,1,function(x){apply(PW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_two)*(x+y+as.vector(AAW)))})}))

         AA_no_PW_3 <- t(apply(PW_all_no,1,function(x){apply(PW_one_no,1,function(y){ sum(All_add(x,y))})}))
         AA_yes_PW_3 <- t(apply(PW_all_yes,1,function(x){apply(PW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_two)*(x+y+as.vector(AAW)))})}))

         AA_no_PW_3t <- t(AA_no_PW_3)
         AA_yes_PW_3t <- t(AA_yes_PW_3)

         AA_no_PW_4 <- t(apply(PW_all_no,1,function(x){apply(PW_all_no,1,function(y){ sum(All_add(x,y))/2})}))
         AA_yes_PW_4 <- t(apply(PW_all_yes,1,function(x){apply(PW_all_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_two)*(x+y+as.vector(AAW))) /2})}))


         CC_no_QW_2 <- t(apply(QW_two_no,1,function(x){apply(QW_one_no,1,function(y){ sum(All_add(x,y))})}))
         CC_yes_QW_2 <- t(apply(QW_two_yes,1,function(x){apply(QW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_two)*(x+y+as.vector(CCW)))})}))

         CC_no_QW_3 <- t(apply(QW_all_no,1,function(x){apply(QW_one_no,1,function(y){ sum(All_add(x,y))})}))
         CC_yes_QW_3 <- t(apply(QW_all_yes,1,function(x){apply(QW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_two)*(x+y+as.vector(CCW)))})}))

         CC_no_QW_4 <- t(apply(QW_all_no,1,function(x){apply(QW_all_no,1,function(y){ sum(All_add(x,y))/2 })}))
         CC_yes_QW_4 <- t(apply(QW_all_yes,1,function(x){apply(QW_all_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_two)*(x+y+as.vector(CCW)))/2 })}))





         M122_1_1 <- mat_tri((A_2BW * CC_no_Q_2 + CC_no_QW_2 * Ubb) * BB_no)
         M122_1_2 <- mat_tri((A_2BW * CC_no_Q_2 + CC_no_QW_2 * Ubb) * BB_two + Ubb * CC_no_Q_2 * BBW)
         M122_1_3 <- mat_tri((A_2BW * CC_yes_Q_2 + CC_yes_QW_2 * Ubb) * BB_no)
         M122_1_4 <- mat_tri((A_2BW * CC_yes_Q_2 + CC_yes_QW_2 * Ubb) * BB_two + Ubb * CC_yes_Q_2 * BBW)

         M122_2_1 <- sum((A_2BW * CC_no_Q_3 + CC_no_QW_3 * Ubb) * BB_no)
         M122_2_2 <- sum((A_2BW * CC_no_Q_3 + CC_no_QW_3 * Ubb) * BB_two + Ubb * CC_no_Q_3 * BBW)
         M122_2_3 <- sum((A_2BW * CC_yes_Q_3 + CC_yes_QW_3 * Ubb) * BB_no)
         M122_2_4 <- sum((A_2BW * CC_yes_Q_3 + CC_yes_QW_3 * Ubb) * BB_two + Ubb * CC_yes_Q_3 * BBW)

         M122_3_1 <- mat_tri((A_2BW * CC_no_Q_4 + CC_no_QW_4 * Ubb) * BB_no)
         M122_3_2 <- mat_tri((A_2BW * CC_no_Q_4 + CC_no_QW_4 * Ubb) * BB_two + Ubb * CC_no_Q_4 * BBW)
         M122_3_3 <- mat_tri((A_2BW * CC_yes_Q_4 + CC_yes_QW_4 * Ubb) * BB_no)
         M122_3_4 <- mat_tri((A_2BW * CC_yes_Q_4 + CC_yes_QW_4 * Ubb) * BB_two + Ubb * CC_yes_Q_4 * BBW)


         M221_1_1 <- mat_tri((AA_no_PW_2 * Vbb + C_2BW * AA_no_P_2) * BB_no)
         M221_1_2 <- mat_tri((AA_yes_PW_2 * Vbb + C_2BW * AA_yes_P_2) * BB_no)
         M221_1_3 <- mat_tri((AA_no_PW_2 * Vbb + C_2BW * AA_no_P_2) * BB_two + AA_no_P_2 * Vbb * BBW)
         M221_1_4 <- mat_tri((AA_yes_PW_2 * Vbb + C_2BW * AA_yes_P_2) * BB_two + AA_yes_P_2 * Vbb * BBW)

         M221_2_1 <- sum((AA_no_PW_3 * Vbb + C_2BW * AA_no_P_3) * BB_no)
         M221_2_2 <- sum((AA_yes_PW_3 * Vbb + C_2BW * AA_yes_P_3) * BB_no)
         M221_2_3 <- sum((AA_no_PW_3 * Vbb + C_2BW * AA_no_P_3) * BB_two + AA_no_P_3 * Vbb * BBW)
         M221_2_4 <- sum((AA_yes_PW_3 * Vbb + C_2BW * AA_yes_P_3) * BB_two + AA_yes_P_3 * Vbb * BBW)

         M221_3_1 <- mat_tri((AA_no_PW_4 * Vbb + C_2BW * AA_no_P_4) * BB_no)
         M221_3_2 <- mat_tri((AA_yes_PW_4 * Vbb + C_2BW * AA_yes_P_4) * BB_no)
         M221_3_3 <- mat_tri((AA_no_PW_4 * Vbb + C_2BW * AA_no_P_4) * BB_two + AA_no_P_4 * Vbb * BBW)
         M221_3_4 <- mat_tri((AA_yes_PW_4 * Vbb + C_2BW * AA_yes_P_4) * BB_two + AA_yes_P_4 * Vbb * BBW)



         M222_1_1 <- sum((AA_no_PW_3 * CC_no_Q_2 + CC_no_QW_2 *AA_no_P_3 ) * BB_no )
         M222_1_2 <- sum((AA_yes_PW_3 * CC_no_Q_2 + CC_no_QW_2 *AA_yes_P_3 ) * BB_no)
         M222_1_3 <- sum((AA_no_PW_3 * CC_no_Q_2 + CC_no_QW_2 *AA_no_P_3) * BB_two + AA_no_P_3 * CC_no_Q_2 * BBW )
         M222_1_4 <- sum((AA_no_PW_3 * CC_yes_Q_2 + CC_yes_QW_2 *AA_no_P_3 ) * BB_no)
         M222_1_5 <- sum((AA_yes_PW_3 * CC_no_Q_2 + CC_no_QW_2 *AA_yes_P_3) * BB_two + AA_yes_P_3 * CC_no_Q_2 * BBW )
         M222_1_6 <- sum((AA_yes_PW_3 * CC_yes_Q_2 + CC_yes_QW_2 *AA_yes_P_3 ) * BB_no)
         M222_1_7 <- sum((AA_no_PW_3 * CC_yes_Q_2 + CC_yes_QW_2 *AA_no_P_3) * BB_two + AA_no_P_3 * CC_yes_Q_2 * BBW )
         M222_1_8 <- sum((AA_yes_PW_3 * CC_yes_Q_2 + CC_yes_QW_2 *AA_yes_P_3) * BB_two + AA_yes_P_3 * CC_yes_Q_2 * BBW )

         M222_2_1 <- mat_tri((AA_no_PW_4 * CC_no_Q_2 + CC_no_QW_2 *AA_no_P_4 ) * BB_no )
         M222_2_2 <- mat_tri((AA_yes_PW_4 * CC_no_Q_2 + CC_no_QW_2 *AA_yes_P_4 ) * BB_no)
         M222_2_3 <- mat_tri((AA_no_PW_4 * CC_no_Q_2 + CC_no_QW_2 *AA_no_P_4) * BB_two + AA_no_P_4 * CC_no_Q_2 * BBW )
         M222_2_4 <- mat_tri((AA_no_PW_4 * CC_yes_Q_2 + CC_yes_QW_2 *AA_no_P_4 ) * BB_no)
         M222_2_5 <- mat_tri((AA_yes_PW_4 * CC_no_Q_2 + CC_no_QW_2 *AA_yes_P_4) * BB_two + AA_yes_P_4 * CC_no_Q_2 * BBW )
         M222_2_6 <- mat_tri((AA_yes_PW_4 * CC_yes_Q_2 + CC_yes_QW_2 *AA_yes_P_4 ) * BB_no)
         M222_2_7 <- mat_tri((AA_no_PW_4 * CC_yes_Q_2 + CC_yes_QW_2 *AA_no_P_4) * BB_two + AA_no_P_4 * CC_yes_Q_2 * BBW )
         M222_2_8 <- mat_tri((AA_yes_PW_4 * CC_yes_Q_2 + CC_yes_QW_2 *AA_yes_P_4) * BB_two + AA_yes_P_4 * CC_yes_Q_2 * BBW )


         M222_3_1 <- sum((AA_no_PW_2 * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_2 ) * BB_no )
         M222_3_2 <- sum((AA_yes_PW_2 * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_2 ) * BB_no)
         M222_3_3 <- sum((AA_no_PW_2 * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_2) * BB_two + AA_no_P_2 * CC_no_Q_3 * BBW )
         M222_3_4 <- sum((AA_no_PW_2 * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_2 ) * BB_no)
         M222_3_5 <- sum((AA_yes_PW_2 * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_2) * BB_two + AA_yes_P_2 * CC_no_Q_3 * BBW )
         M222_3_6 <- sum((AA_yes_PW_2 * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_2 ) * BB_no)
         M222_3_7 <- sum((AA_no_PW_2 * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_2) * BB_two + AA_no_P_2 * CC_yes_Q_3 * BBW )
         M222_3_8 <- sum((AA_yes_PW_2 * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_2) * BB_two + AA_yes_P_2 * CC_yes_Q_3 * BBW )


         M222_4_1 <- sum((AA_no_PW_3 * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_3 ) * BB_no )
         M222_4_2 <- sum((AA_yes_PW_3 * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_3 ) * BB_no)
         M222_4_3 <- sum((AA_no_PW_3 * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_3) * BB_two + AA_no_P_3 * CC_no_Q_3 * BBW )
         M222_4_4 <- sum((AA_no_PW_3 * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_3 ) * BB_no)
         M222_4_5 <- sum((AA_yes_PW_3 * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_3) * BB_two + AA_yes_P_3 * CC_no_Q_3 * BBW )
         M222_4_6 <- sum((AA_yes_PW_3 * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_3 ) * BB_no)
         M222_4_7 <- sum((AA_no_PW_3 * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_3) * BB_two + AA_no_P_3 * CC_yes_Q_3 * BBW )
         M222_4_8 <- sum((AA_yes_PW_3 * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_3) * BB_two + AA_yes_P_3 * CC_yes_Q_3 * BBW )

         M222_5_1 <- sum((AA_no_PW_4 * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_4 ) * BB_no )
         M222_5_2 <- sum((AA_yes_PW_4 * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_4 ) * BB_no)
         M222_5_3 <- sum((AA_no_PW_4 * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_4) * BB_two + AA_no_P_4 * CC_no_Q_3 * BBW )
         M222_5_4 <- sum((AA_no_PW_4 * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_4 ) * BB_no)
         M222_5_5 <- sum((AA_yes_PW_4 * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_4) * BB_two + AA_yes_P_4 * CC_no_Q_3 * BBW )
         M222_5_6 <- sum((AA_yes_PW_4 * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_4 ) * BB_no)
         M222_5_7 <- sum((AA_no_PW_4 * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_4) * BB_two + AA_no_P_4 * CC_yes_Q_3 * BBW )
         M222_5_8 <- sum((AA_yes_PW_4 * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_4) * BB_two + AA_yes_P_4 * CC_yes_Q_3 * BBW )

         M222_6_1 <- mat_tri((AA_no_PW_2 * CC_no_Q_4 + CC_no_QW_4 *AA_no_P_2 ) * BB_no )
         M222_6_2 <- mat_tri((AA_yes_PW_2 * CC_no_Q_4 + CC_no_QW_4 *AA_yes_P_2 ) * BB_no)
         M222_6_3 <- mat_tri((AA_no_PW_2 * CC_no_Q_4 + CC_no_QW_4 *AA_no_P_2) * BB_two + AA_no_P_2 * CC_no_Q_4 * BBW )
         M222_6_4 <- mat_tri((AA_no_PW_2 * CC_yes_Q_4 + CC_yes_QW_4 *AA_no_P_2 ) * BB_no)
         M222_6_5 <- mat_tri((AA_yes_PW_2 * CC_no_Q_4 + CC_no_QW_4 *AA_yes_P_2) * BB_two + AA_yes_P_2 * CC_no_Q_4 * BBW )
         M222_6_6 <- mat_tri((AA_yes_PW_2 * CC_yes_Q_4 + CC_yes_QW_4 *AA_yes_P_2 ) * BB_no)
         M222_6_7 <- mat_tri((AA_no_PW_2 * CC_yes_Q_4 + CC_yes_QW_4 *AA_no_P_2) * BB_two + AA_no_P_2 * CC_yes_Q_4 * BBW )
         M222_6_8 <- mat_tri((AA_yes_PW_2 * CC_yes_Q_4 + CC_yes_QW_4 *AA_yes_P_2) * BB_two + AA_yes_P_2 * CC_yes_Q_4 * BBW )

         M222_7_1 <- sum((AA_no_PW_3 * CC_no_Q_4 + CC_no_QW_4 *AA_no_P_3 ) * BB_no )
         M222_7_2 <- sum((AA_yes_PW_3 * CC_no_Q_4 + CC_no_QW_4 *AA_yes_P_3 ) * BB_no)
         M222_7_3 <- sum((AA_no_PW_3 * CC_no_Q_4 + CC_no_QW_4 *AA_no_P_3) * BB_two + AA_no_P_3 * CC_no_Q_4 * BBW )
         M222_7_4 <- sum((AA_no_PW_3 * CC_yes_Q_4 + CC_yes_QW_4 *AA_no_P_3 ) * BB_no)
         M222_7_5 <- sum((AA_yes_PW_3 * CC_no_Q_4 + CC_no_QW_4 *AA_yes_P_3) * BB_two + AA_yes_P_3 * CC_no_Q_4 * BBW )
         M222_7_6 <- sum((AA_yes_PW_3 * CC_yes_Q_4 + CC_yes_QW_4 *AA_yes_P_3 ) * BB_no)
         M222_7_7 <- sum((AA_no_PW_3 * CC_yes_Q_4 + CC_yes_QW_4 *AA_no_P_3) * BB_two + AA_no_P_3 * CC_yes_Q_4 * BBW )
         M222_7_8 <- sum((AA_yes_PW_3 * CC_yes_Q_4 + CC_yes_QW_4 *AA_yes_P_3) * BB_two + AA_yes_P_3 * CC_yes_Q_4 * BBW )

         M222_8_1 <- mat_tri((AA_no_PW_4 * CC_no_Q_4 + CC_no_QW_4 *AA_no_P_4 ) * BB_no )
         M222_8_2 <- mat_tri((AA_yes_PW_4 * CC_no_Q_4 + CC_no_QW_4 *AA_yes_P_4 ) * BB_no)
         M222_8_3 <- mat_tri((AA_no_PW_4 * CC_no_Q_4 + CC_no_QW_4 *AA_no_P_4) * BB_two + AA_no_P_4 * CC_no_Q_4 * BBW )
         M222_8_4 <- mat_tri((AA_no_PW_4 * CC_yes_Q_4 + CC_yes_QW_4 *AA_no_P_4 ) * BB_no)
         M222_8_5 <- mat_tri((AA_yes_PW_4 * CC_no_Q_4 + CC_no_QW_4 *AA_yes_P_4) * BB_two + AA_yes_P_4 * CC_no_Q_4 * BBW )
         M222_8_6 <- mat_tri((AA_yes_PW_4 * CC_yes_Q_4 + CC_yes_QW_4 *AA_yes_P_4 ) * BB_no)
         M222_8_7 <- mat_tri((AA_no_PW_4 * CC_yes_Q_4 + CC_yes_QW_4 *AA_no_P_4) * BB_two + AA_no_P_4 * CC_yes_Q_4 * BBW )
         M222_8_8 <- mat_tri((AA_yes_PW_4 * CC_yes_Q_4 + CC_yes_QW_4 *AA_yes_P_4) * BB_two + AA_yes_P_4 * CC_yes_Q_4 * BBW )

         M222_9_1 <- sum((AA_no_PW_3t * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_3t ) * BB_no )
         M222_9_2 <- sum((AA_yes_PW_3t * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_3t ) * BB_no)
         M222_9_3 <- sum((AA_no_PW_3t * CC_no_Q_3 + CC_no_QW_3 *AA_no_P_3t) * BB_two + AA_no_P_3t * CC_no_Q_3 * BBW )
         M222_9_4 <- sum((AA_no_PW_3t * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_3t ) * BB_no)
         M222_9_5 <- sum((AA_yes_PW_3t * CC_no_Q_3 + CC_no_QW_3 *AA_yes_P_3t) * BB_two + AA_yes_P_3t * CC_no_Q_3 * BBW )
         M222_9_6 <- sum((AA_yes_PW_3t * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_3t ) * BB_no)
         M222_9_7 <- sum((AA_no_PW_3t * CC_yes_Q_3 + CC_yes_QW_3 *AA_no_P_3t) * BB_two + AA_no_P_3t * CC_yes_Q_3 * BBW )
         M222_9_8 <- sum((AA_yes_PW_3t * CC_yes_Q_3 + CC_yes_QW_3 *AA_yes_P_3t) * BB_two + AA_yes_P_3t * CC_yes_Q_3 * BBW )

         motif_weighted <- c(M111,M112_1,M112_2,M211_1,M211_2,M212_1,M212_2,M212_3,M212_4,M121_1,M121_2,M122_1_1,M122_1_2,M122_1_3,M122_1_4,
                             M122_2_1,M122_2_2,M122_2_3,M122_2_4,M122_3_1,M122_3_2,M122_3_3,M122_3_4,M221_1_1,M221_1_2,M221_1_3,M221_1_4,
                             M221_2_1,M221_2_2,M221_2_3,M221_2_4,M221_3_1,M221_3_2,M221_3_3,M221_3_4,M222_1_1,M222_1_2,M222_1_3,M222_1_4,
                             M222_1_5,M222_1_6,M222_1_7,M222_1_8,M222_2_1,M222_2_2,M222_2_3,M222_2_4,M222_2_5,M222_2_6,M222_2_7,M222_2_8,
                             M222_3_1,M222_3_2,M222_3_3,M222_3_4,M222_3_5,M222_3_6,M222_3_7,M222_3_8,
                             M222_4_1,M222_4_2,M222_4_3,M222_4_4,M222_4_5,M222_4_6,M222_4_7,M222_4_8,
                             M222_5_1,M222_5_2,M222_5_3,M222_5_4,M222_5_5,M222_5_6,M222_5_7,M222_5_8,
                             M222_6_1,M222_6_2,M222_6_3,M222_6_4,M222_6_5,M222_6_6,M222_6_7,M222_6_8,
                             M222_7_1,M222_7_2,M222_7_3,M222_7_4,M222_7_5,M222_7_6,M222_7_7,M222_7_8,
                             M222_8_1,M222_8_2,M222_8_3,M222_8_4,M222_8_5,M222_8_6,M222_8_7,M222_8_8,
                             M222_9_1,M222_9_2,M222_9_3,M222_9_4,M222_9_5,M222_9_6,M222_9_7,M222_9_8)

         intra_guild_edge <- c(2, 3, 4, 3, 4, 4, 5, 5, 6, 4, 5, 4, 5, 5, 6, 5, 6, 6, 7, 6, 7, 7, 8, 4, 5, 5, 6, 5, 6, 6, 7, 6, 7, 7,
         8, 5, 6, 6, 6, 7, 7, 7, 8, 6, 7, 7, 7, 8, 8, 8, 9, 5, 6, 6, 6, 7, 7, 7, 8, 6, 7, 7, 7, 8, 8, 8, 9, 7, 8, 8, 8, 9, 9, 9,
         10, 6, 7, 7, 7, 8, 8, 8, 9, 7, 8, 8, 8, 9, 9, 9, 10, 8, 9, 9, 9, 10, 10, 10, 11, 6, 7, 7, 7, 8, 8, 8, 9)

         mean_weighted <- motif_weighted/(intra_guild_motif * intra_guild_edge)

         mean_weighted <- replace(mean_weighted,which(is.nan(mean_weighted)),NA)


         return(data.frame(motif_name=motif_names,count=intra_guild_motif,mean_weight=mean_weighted ))
      }
      else{
         return(data.frame(motif_name=motif_names,count=intra_guild_motif))
      }
   }
