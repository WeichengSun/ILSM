#' Count the roles of connector nodes defined by interconnection motifs with intra-guild interactions
#'
#' Count the roles of connector nodes defined by interconnection motifs in a tripartite network with intra-guild interactions.
#'
#' @param mat A square block interaction matrix representing a tripartite network including intra-guild and inter-guild interactions. See details.
#' @param guilds A character vector matching rows of \code{mat} to indicate the guilds using ('a','b' and 'c'). See details.
#' @param weighted Logical. Default to FALSE. If TRUE, a weighted measure is provided. See details.
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
#' \code{guilds} should be a vector of the same length as the row of \code{mat} like c("a","a"..."b","b"..."c","c"..)
#'
#' \strong{Interconnection motifs in tripartite networks with intra-guild interactions}
#' <br>An interconnection motif is defined to comprise three sets of connected nodes: the connector nodes (belonging to b-nodes), the nodes in one subnetwork (belonging to a-nodes in the P subnetwork), and the nodes in the other subnetwork (belonging to c-nodes in the Q subnetwork). Each guild has two nodes at most, resulting in a total of 107 distinct motif forms.
#' The algorithm for counting interconnection motifs is designed by extending the fast approach from Simmons et al.(2019). Currently, only the roles in first 35 motifs are provided (see ). For interconnection motifs in tripartite networks without intra-guild interactions, please see **icmotif_count** and **icmotif_role**.
#'
#' \strong{Weighted networks}
#' <br>For weighted tripartite networks, the mean weight of the motif occurrence (i.e., a motif occurrence isomorphic to a particular motif form) is provided for a given node with a given role, following Mora et al. (2018) and Simmons et al. (2019).
#'
#' @return
#' For binary networks, return a matrix with elements representing the number of times each connector node plays for each unique role within interconnection motifs; for weighted networks, the matrix element represents the mean weight of the motif occurrences where the node exist.
#'
#'
#' @references
#' García-Callejas, D., Godoy, O., Buche, L., Hurtado, M., Lanuza, J.B., Allen-Perkins, A. et al. (2023) Non-random interactions within and across guilds shape the potential to coexist in multi-trophic ecological communities. Ecology Letters, 26, 831–842.
#'
#' Pilosof, S., Porter, M. A., Pascual, M., & KC)fi, S. (2017). The multilayer nature of ecological networks. Nature Ecology & Evolution, 1(4), 0101.
#'
#' Mora, B.B., Cirtwill, A.R. and Stouffer, D.B. (2018). pymfinder: a tool for the motif analysis of binary and quantitative complex networks. bioRxiv, 364703.
#'
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#' @examples
#'
#' ## A toy tripartite network with intra-guild negative interactions, inter-guild mutualistic interactions and inter-guild antagonistic interactions.
#' set.seed(12)
#' ##4 a-nodes,5 b-nodes, and 3 c-nodes
#'
#' ##intra-guild interaction matrices
#' mat_aa<-matrix(runif(16,-0.8,-0.2),4,4)
#' mat_aa <- mat_aa + t(mat_aa); diag(mat_aa) <- 0
#' mat_bb<-matrix(runif(25,-0.8,-0.2),5,5)
#' mat_bb <- mat_bb + t(mat_bb); diag(mat_bb) <- 0
#' mat_cc<-matrix(runif(9,-0.8,-0.2),3,3)
#' mat_cc <- mat_cc + t(mat_cc); diag(mat_cc) <- 0
#'
#' ##inter-guild interaction matrices between a- and b-nodes.
#' mat_ab<-mat_ba<-matrix(sample(c(rep(0,8),runif(12,0,0.5))),4,5,byrow=T)# interaction probability = 12/20
#' mat_ba<-t(mat_ba)
#'
#' ##inter-guild interaction matrices between b- and c-nodes.
#' mat_cb<-mat_bc<-matrix(sample(c(rep(0,8),runif(7,0,0.5))),3,5,byrow=T)# interaction probability = 7/15
#' mat_bc<--t(mat_bc)
#' mat<-rbind(cbind(mat_aa,mat_ab,matrix(0,4,3)),cbind(mat_ba,mat_bb,mat_bc),cbind(matrix(0,3,4),mat_cb,mat_cc))
#'
#' ##set the node names
#' rownames(mat)<-c(paste0("a",1:4),paste0("b",1:5),paste0("c",1:3));colnames(mat)<-c(paste0("a",1:4),paste0("b",1:5),paste0("c",1:3))
#'
#' myguilds=c(rep("a",4),rep("b",5),rep("c",3))
#' ig_icmotif_role(mat,guilds=myguilds)
#' ig_icmotif_role(mat,guilds=myguilds,T)
#'
#'
#'
ig_icmotif_role<-function(mat, guilds,weighted =FALSE) {
      if (inherits(mat, c("matrix", "data.frame")) == T || ncol(mat)==nrow(mat)){
         if(length(guilds)!= nrow(mat))
            stop("guilds should be a vector of the same length as the row of mat!")
         else{
            if(sum(unique(guilds)%in%c("a","b","c"))!=3)
               stop("guilds should be composed of these three characters (a,b,c)!")
            else{
               if(is.null(dimnames(mat)) || sum(rownames(mat)!=colnames(mat))!=0 )
                  stop("Make sure this matrix either have the same row and column names [No NA!!!]")
               else{
                  mat <- abs(mat)
                  subnet_mat1 <- mat[guilds == "a", guilds == "b"]
                  subnet_mat2 <- mat[guilds == "b", guilds == "c"]

                  AA <- mat[guilds == "a", guilds == "a"]

                  BB <- mat[guilds == "b", guilds == "b"]

                  CC <- mat[guilds == "c", guilds == "c"]
               }
            }
         }
      }
      else
         stop("mat should be a square block interaction matrix!")

      dimnames(subnet_mat1) <- NULL
      dimnames(subnet_mat2) <- NULL
      P <- subnet_mat1
      Q <- subnet_mat2

      PW <- P; QW <- Q; AAW <- AA; BBW <- BB; CCW <- CC
      AA <- (AA!=0)*1; BB <- (BB!=0)*1; CC <- (CC!=0)*1
      P <- (P!=0)*1; Q <- (Q!=0)*1

      AA_no <- AA == 0
      diag(AA) <- 0
      diag(AA_no) <- 0
      AA_tow <- (AA!=0)*1

      BB_no <- BB == 0
      diag(BB) <- 0
      diag(BB_no) <- 0
      BB_two <- (BB!=0)*1


      CC_no <- CC == 0
      diag(CC) <- 0
      diag(CC_no) <- 0
      CC_tow <- (CC!=0)*1


      Ob <- colSums(P); Rb <- rowSums(Q)
      Ubb <- t(P) %*% P
      Vbb <- Q %*% t(Q)
      UVbb <- Ubb * Vbb
      B2C_NO <- rowSums((Q%*%CC_no) *Q) /2
      B2C_YES <- rowSums((Q%*%CC_tow) *Q) /2
      B2A_NO <- rowSums((t(P)%*%AA_no) *t(P)) /2
      B2A_YES <- rowSums((t(P)%*%AA_tow) *t(P)) /2

      BB_role <- NULL

      M111 <- Ob * Rb
      BB_role <- cbind(BB_role,M111)

      M112_1 <- Ob * B2C_NO
      BB_role <- cbind(BB_role,M112_1)

      M112_2 <- Ob * B2C_YES
      BB_role <- cbind(BB_role,M112_2)

      M211_1 <- B2A_NO * Rb
      BB_role <- cbind(BB_role,M211_1)

      M211_2 <- B2A_YES * Rb
      BB_role <- cbind(BB_role,M211_2)

      M212_1 <- B2A_NO * B2C_NO
      BB_role <- cbind(BB_role,M212_1)

      M212_2 <- B2A_YES * B2C_NO
      BB_role <- cbind(BB_role,M212_2)

      M212_3 <- B2A_NO * B2C_YE
      BB_role <- cbind(BB_role,M212_3)

      M212_4 <- B2A_YES * B2C_YES
      BB_role <- cbind(BB_role,M212_4)

      # mat_tri <- function(mat) {sum(mat[upper.tri(mat)])}

      M121_1 <- rowSums(UVbb * BB_no)
      BB_role <- cbind(BB_role,M121_1)

      M121_2 <- rowSums(UVbb * BB_two)
      BB_role <- cbind(BB_role,M121_2)


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



      M122_1_1 <- rowSums(Ubb * BB_no * CC_no_Q_2)
      BB_role <- cbind(BB_role,M122_1_1)

      M122_1_2 <- rowSums(Ubb * BB * CC_no_Q_2)
      BB_role <- cbind(BB_role,M122_1_2)

      M122_1_3 <- rowSums(Ubb * BB_no * CC_yes_Q_2)
      BB_role <- cbind(BB_role,M122_1_3)

      M122_1_4 <- rowSums(Ubb * BB * CC_yes_Q_2)
      BB_role <- cbind(BB_role,M122_1_4)

      M122_2_1 <- Ubb * BB_no * CC_no_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_1), colSums(M122_2_1))

      M122_2_2 <- Ubb * BB * CC_no_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_2), colSums(M122_2_2))

      M122_2_3 <- Ubb * BB_no * CC_yes_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_3), colSums(M122_2_3))

      M122_2_4 <- Ubb * BB * CC_yes_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_4), colSums(M122_2_4))

      M122_3_1 <- rowSums(Ubb * BB_no * CC_no_Q_4)
      BB_role <- cbind(BB_role,M122_3_1)

      M122_3_2 <- rowSums(Ubb * BB * CC_no_Q_4)
      BB_role <- cbind(BB_role,M122_3_2)

      M122_3_3 <- rowSums(Ubb * BB_no * CC_yes_Q_4)
      BB_role <- cbind(BB_role,M122_3_3)

      M122_3_4 <- rowSums(Ubb * BB * CC_yes_Q_4)
      BB_role <- cbind(BB_role,M122_3_4)

      M221_1_1 <- rowSums(AA_no_P_2 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_1_1)

      M221_1_2 <- rowSums(AA_yes_P_2 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_1_2)

      M221_1_3 <- rowSums(AA_no_P_2 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_1_3)

      M221_1_4 <- rowSums(AA_yes_P_2 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_1_4)


      M221_2_1 <- sum(AA_no_P_3 * BB_no * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_1), colSums(M221_2_1))

      M221_2_2 <- sum(AA_yes_P_3 * BB_no * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_2), colSums(M221_2_2))

      M221_2_3 <- sum(AA_no_P_3 * BB * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_3), colSums(M221_2_3))

      M221_2_4 <- sum(AA_yes_P_3 * BB * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_4), colSums(M221_2_4))

      M221_3_1 <- rowSums(AA_no_P_4 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_3_1)

      M221_3_2 <- rowSums(AA_yes_P_4 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_3_2)

      M221_3_3 <- rowSums(AA_no_P_4 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_3_3)

      M221_3_4 <- rowSums(AA_yes_P_4 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_3_4)

      rownames(BB_role) <- matrow
      colnames(BB_role) <- paste0("role",1:ncol(BB_role))

      intra_guild_role <- BB_role



      if(weighted){
         All_add <- function(a,b){ return((a!=0)*(b!=0)*(a+b)) }
         One_add <- function(a,b){ return((a!=0)*(b!=0)*(a)) }

         ObW <- colSums(PW)
         RbW <- rowSums(QW)

         B2C_NOW <- apply(QW,1,function(x){ A<- apply(CC_no,2,function(y){ One_add(x,y) });  return(sum(All_add(x,t(A)))/2)})
         B2C_YESW <- apply(QW,1,function(x){ A<- apply(CCW,2,function(y){ All_add(x,y) });  return(sum(All_add(x,t(A)))/2)})

         B2A_NOW <- apply(PW,2,function(x){ A<- apply(AA_no,2,function(y){ One_add(x,y) });  return(sum(All_add(x,t(A)))/2)})
         B2A_YESW <- apply(PW,2,function(x){ A<- apply(AAW,2,function(y){ All_add(x,y) });  return(sum(All_add(x,t(A)))/2)})

         BB_role <- NULL

         M111 <- (ObW * Rb + RbW * Ob)/2
         BB_role <- cbind(BB_role,M111)

         M112_1 <- (ObW * B2C_NO + B2C_NOW * Ob  )/3
         BB_role <- cbind(BB_role,M112_1)

         M112_2 <- (ObW * B2C_YES + B2C_YESW * Ob  )/4
         BB_role <- cbind(BB_role,M112_2)

         M211_1 <- (B2A_NOW * Rb + RbW * B2A_NO )/3
         BB_role <- cbind(BB_role,M211_1)

         M211_2 <- (B2A_YESW * Rb +  RbW * B2A_YES )/4
         BB_role <- cbind(BB_role,M211_2)

         M212_1 <- (B2A_NOW * B2C_NO + B2C_NOW *B2A_NO)/4
         BB_role <- cbind(BB_role,M212_1)

         M212_2 <- (B2A_YESW * B2C_NO + B2C_NOW * B2A_YES)/5
         BB_role <- cbind(BB_role,M212_2)

         M212_3 <- (B2A_NOW * B2C_YES + B2C_YESW * B2A_NO)/5
         BB_role <- cbind(BB_role,M212_3)

         M212_4 <- (B2A_YESW * B2C_YES + B2C_YESW * B2A_YES)/6
         BB_role <- cbind(BB_role,M212_4)


         A_2BW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum(All_add(x,y)) })}))
         C_2BW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum(All_add(x,y)) })}))


         M121_1 <- rowSums( (A_2BW*Vbb + C_2BW*Ubb) * BB_no)/4
         BB_role <- cbind(BB_role,M121_1)

         M121_2 <- rowSums((A_2BW*Vbb + C_2BW*Ubb)* BB_two + UVbb * BBW )/5
         BB_role <- cbind(BB_role,M121_2)


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
         AA_yes_PW_2 <- t(apply(PW_two_yes,1,function(x){apply(PW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_tow)*(x+y+as.vector(AAW)))})}))

         AA_no_PW_3 <- t(apply(PW_all_no,1,function(x){apply(PW_one_no,1,function(y){ sum(All_add(x,y))})}))
         AA_yes_PW_3 <- t(apply(PW_all_yes,1,function(x){apply(PW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_tow)*(x+y+as.vector(AAW)))})}))

         AA_no_PW_3t <- t(AA_no_PW_3)
         AA_yes_PW_3t <- t(AA_yes_PW_3)

         AA_no_PW_4 <- t(apply(PW_all_no,1,function(x){apply(PW_all_no,1,function(y){ sum(All_add(x,y))/2})}))
         AA_yes_PW_4 <- t(apply(PW_all_yes,1,function(x){apply(PW_all_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_tow)*(x+y+as.vector(AAW))) /2})}))


         CC_no_QW_2 <- t(apply(QW_two_no,1,function(x){apply(QW_one_no,1,function(y){ sum(All_add(x,y))})}))
         CC_yes_QW_2 <- t(apply(QW_two_yes,1,function(x){apply(QW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_tow)*(x+y+as.vector(CCW)))})}))

         CC_no_QW_3 <- t(apply(QW_all_no,1,function(x){apply(QW_one_no,1,function(y){ sum(All_add(x,y))})}))
         CC_yes_QW_3 <- t(apply(QW_all_yes,1,function(x){apply(QW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_tow)*(x+y+as.vector(CCW)))})}))

         CC_no_QW_4 <- t(apply(QW_all_no,1,function(x){apply(QW_all_no,1,function(y){ sum(All_add(x,y))/2 })}))
         CC_yes_QW_4 <- t(apply(QW_all_yes,1,function(x){apply(QW_all_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_tow)*(x+y+as.vector(CCW)))/2 })}))





         M122_1_1 <- rowSums((A_2BW * CC_no_Q_2 + CC_no_QW_2 * Ubb) * BB_no)/4
         BB_role <- cbind(BB_role,M122_1_1)

         M122_1_2 <- rowSums((A_2BW * CC_no_Q_2 + CC_no_QW_2 * Ubb) * BB_two + Ubb * CC_no_Q_2 * BBW)/5
         BB_role <- cbind(BB_role,M122_1_2)

         M122_1_3 <- rowSums((A_2BW * CC_yes_Q_2 + CC_yes_QW_2 * Ubb) * BB_no)/5
         BB_role <- cbind(BB_role,M122_1_3)

         M122_1_4 <- rowSums((A_2BW * CC_yes_Q_2 + CC_yes_QW_2 * Ubb) * BB_two + Ubb * CC_yes_Q_2 * BBW)/6
         BB_role <- cbind(BB_role,M122_1_4)

         M122_2_1 <- ((A_2BW * CC_no_Q_3 + CC_no_QW_3 * Ubb) * BB_no)/5
         BB_role <- cbind(BB_role,rowSums(M122_2_1), colSums(M122_2_1))

         M122_2_2 <- ((A_2BW * CC_no_Q_3 + CC_no_QW_3 * Ubb) * BB_two + Ubb * CC_no_Q_3 * BBW)/6
         BB_role <- cbind(BB_role,rowSums(M122_2_2), colSums(M122_2_2))

         M122_2_3 <- ((A_2BW * CC_yes_Q_3 + CC_yes_QW_3 * Ubb) * BB_no)/6
         BB_role <- cbind(BB_role,rowSums(M122_2_3), colSums(M122_2_3))

         M122_2_4 <- ((A_2BW * CC_yes_Q_3 + CC_yes_QW_3 * Ubb) * BB_two + Ubb * CC_yes_Q_3 * BBW)/7
         BB_role <- cbind(BB_role,rowSums(M122_2_4), colSums(M122_2_4))



         M122_3_1 <- rowSums((A_2BW * CC_no_Q_4 + CC_no_QW_4 * Ubb) * BB_no)/6
         BB_role <- cbind(BB_role,M122_3_1)

         M122_3_2 <- rowSums((A_2BW * CC_no_Q_4 + CC_no_QW_4 * Ubb) * BB_two + Ubb * CC_no_Q_4 * BBW)/7
         BB_role <- cbind(BB_role,M122_3_2)

         M122_3_3 <- rowSums((A_2BW * CC_yes_Q_4 + CC_yes_QW_4 * Ubb) * BB_no)/7
         BB_role <- cbind(BB_role,M122_3_3)

         M122_3_4 <- rowSums((A_2BW * CC_yes_Q_4 + CC_yes_QW_4 * Ubb) * BB_two + Ubb * CC_yes_Q_4 * BBW)/8
         BB_role <- cbind(BB_role,M122_3_4)


         M221_1_1 <- rowSums((AA_no_PW_2 * Vbb + C_2BW * AA_no_P_2) * BB_no)/4
         BB_role <- cbind(BB_role,M221_1_1)

         M221_1_2 <- rowSums((AA_yes_PW_2 * Vbb + C_2BW * AA_yes_P_2) * BB_no)/5
         BB_role <- cbind(BB_role,M221_1_2)

         M221_1_3 <- rowSums((AA_no_PW_2 * Vbb + C_2BW * AA_no_P_2) * BB_two + AA_no_P_2 * Vbb * BBW)/5
         BB_role <- cbind(BB_role,M221_1_3)

         M221_1_4 <- rowSums((AA_yes_PW_2 * Vbb + C_2BW * AA_yes_P_2) * BB_two + AA_yes_P_2 * Vbb * BBW)/6
         BB_role <- cbind(BB_role,M221_1_4)



         M221_2_1 <- ((AA_no_PW_3 * Vbb + C_2BW * AA_no_P_3) * BB_no)/5
         BB_role <- cbind(BB_role,rowSums(M221_2_1), colSums(M221_2_1))

         M221_2_2 <- ((AA_yes_PW_3 * Vbb + C_2BW * AA_yes_P_3) * BB_no)/6
         BB_role <- cbind(BB_role,rowSums(M221_2_2), colSums(M221_2_2))

         M221_2_3 <- ((AA_no_PW_3 * Vbb + C_2BW * AA_no_P_3) * BB_two + AA_no_P_3 * Vbb * BBW)/6
         BB_role <- cbind(BB_role,rowSums(M221_2_3), colSums(M221_2_3))

         M221_2_4 <- ((AA_yes_PW_3 * Vbb + C_2BW * AA_yes_P_3) * BB_two + AA_yes_P_3 * Vbb * BBW)/7
         BB_role <- cbind(BB_role,rowSums(M221_2_4), colSums(M221_2_4))


         M221_3_1 <- rowSums((AA_no_PW_4 * Vbb + C_2BW * AA_no_P_4) * BB_no)/6
         BB_role <- cbind(BB_role,M221_3_1)

         M221_3_2 <- rowSums((AA_yes_PW_4 * Vbb + C_2BW * AA_yes_P_4) * BB_no)/7
         BB_role <- cbind(BB_role,M221_3_2)

         M221_3_3 <- rowSums((AA_no_PW_4 * Vbb + C_2BW * AA_no_P_4) * BB_two + AA_no_P_4 * Vbb * BBW)/7
         BB_role <- cbind(BB_role,M221_3_3)

         M221_3_4 <- rowSums((AA_yes_PW_4 * Vbb + C_2BW * AA_yes_P_4) * BB_two + AA_yes_P_4 * Vbb * BBW)/8
         BB_role <- cbind(BB_role,M221_3_4)

         rownames(BB_role) <- matrow
         colnames(BB_role) <- paste0("weighted role",1:ncol(BB_role))
         intra_guild_mean_role <- BB_role

         intra_guild_mean_role <- intra_guild_mean_role/intra_guild_role

         intra_guild_mean_role <- replace(intra_guild_mean_role, which(is.nan(intra_guild_mean_role)),NA)
         return(intra_guild_mean_role)
      }
      else
         return(intra_guild_role)
   }
