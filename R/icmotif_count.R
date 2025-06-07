#' Count interconnection motifs
#'
#' Counting the frequencies of interconnection motifs for a tripartite interaction network.
#'
#' @param network.or.subnet_mat1 An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.
#' @param subnet_mat2 A matrix representing one subnetwork.
#' @param weighted Logical. Default to FALSE. If TRUE, the arithmetic mean of the subgraph weights is provided for each motif. See details

#' @import igraph
#'
#' @export
#'
#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#' <br>An interconnection motif is defined to comprise three sets of connected nodes: the connector nodes (belonging to b-nodes), the nodes in one subnetwork (belonging to a-nodes in the P subnetwork), and the nodes in the other subnetwork (belonging to c-nodes in the Q subnetwork). Each motif has maximumly 6 nodes, resulting in a total of 48 distinct motif forms.
#' The algorithm for counting interconnection motifs is designed by extending the fast approach from Simmons et al.(2019), which uses mathematical operations directly on the bi-adjacency matrix. For interconnection motifs in tripartite networks with intra-guild interactions, please see **ig_icmotif_count** and **ig_icmotif_role**.
#'
#' Two types of inputs \code{network.or.subnet_mat1} can be processed:
#' \itemize{
#' \item An "igraph" object with node attribute 'level' (0 for a-nodes, 1 for b-nodes,2 for c-nodes). If the input is a weighted network, the edge should have a 'weight' attribute.
#' \item Or a matrix representing subnetwork P, and must be input with \code{subnet_mat2} representing subnetwork Q.
#' }
#'
#' If the inputs are two matrices, please make sure the rows of
#'  \code{network.or.subnet_mat1} and \code{subnet_mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
#'  names to produce connector nodes. Otherwise, row numbers are assigned to row names and matched. Within the two matrices (P and Q), columns represents a-group nodes and c-group nodes respectively.
#'  Elements in matrices are non-zero values if two nodes are linked with or without weights, and 0 otherwise.
#'
#' \strong{Weighted networks}
#' <br>For weighted tripartite networks, the mean weight of a given motif is provided by averaging the weights of all the cases of a particular motif. The weight of a motif case is the arithmetic mean of the weights of its links, following Mora et al. (2018) and Simmons et al. (2019).
#'

#'
#' @return
#'  Return a data.fame of the frequencies (and mean weight) of 48 interconnection motifs. See 'Multi_motif' for the forms.
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
#' Net <- build_toy_net(11,15,16,0.2)
#' icmotif_count(Net)
#'
#' ## empirical network
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' icmotif_count(Net)
#' set.seed(13)
#' library(igraph)
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
icmotif_count <- function(network.or.subnet_mat1, subnet_mat2=NULL, weighted =FALSE){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network<-adjust_net(network.or.subnet_mat1)
      PQ <- as.matrix(network[])
      dimnames(PQ)<-NULL
      P <- PQ[(V(network)$level)==0,(V(network)$level)==1]
      Q <- PQ[(V(network)$level)==1,(V(network)$level)==2]
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
      mat1<-network.or.subnet_mat1
      mat2<-subnet_mat2
      if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
         rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
         rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
         matrow<-unique(c(rownames(mat1),rownames(mat2)))
      }
     #if(nrow(mat1)!=nrow(mat2))
     #   message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
      if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
         matrow<-unique(c(rownames(mat1),rownames(mat2)))
      else
         stop("Make sure matrices either have no row names or have full row names. No NA!")
      mat_1<-matrix(0,length(matrow),ncol(mat1))
      rownames(mat_1)<-matrow
      mat_1[rownames(mat1),]<-mat1
      mat_2<-matrix(0,length(matrow),ncol(mat2))
      rownames(mat_2)<-matrow
      mat_2[rownames(mat2),]<-mat2

      P <- t(mat_1); Q <- mat_2
      dimnames(P)<-NULL; dimnames(Q)<-NULL

      logi<-(apply(P,2,sum)*apply(Q,1,sum))!=0
      P <- P[,logi]
      Q <- Q[logi,]
   }
   else
      stop("Please check the data type of network.or.subnet_mat1 and other parameters.")

   L1<-ncol(P)
   if(L1 < 4L)
      stop("Please input a network with the number of 'connector node' >=4.")

   PW <- P;  QW <- Q

   P <- (P!=0)*1; Q <- (Q!=0)*1
   Ob <- colSums(P)
   Rb <- rowSums(Q)
   Ubb <- t(P)%*%P
   Xbb <- t(P)%*%(1-P)
   XT <- t(Xbb)
   Vbb <- Q%*%t(Q)
   Ybb <- Q%*%(1-t(Q))
   YT <- t(Ybb)
   UVbb <- Ubb*Vbb
   Two<-function(a){ return(a*(a-1)/2) }
   Three<-function(a){ return(a*(a-1)*(a-2)/6)}
   Four<-function(a){ return(a*(a-1)*(a-2)*(a-3)/24)}

   PQ <- P%*%Q

   M111 <- sum(Ob*Rb) #=sum(PQ)

   M112 <- sum(Ob*Two(Rb))

   M113 <- sum(Ob*Three(Rb))

   M114 <- sum(Ob*Four(Rb))

   M211 <- sum(Two(Ob)*Rb)

   M212 <- sum(Two(Ob)*Two(Rb))

   M213 <- sum(Two(Ob)*Three(Rb))

   M311 <- sum(Three(Ob)*Rb)

   M312 <- sum(Three(Ob)*Two(Rb))

   M411 <- sum(Four(Ob)*Rb)

   M121 <- (sum(UVbb)-M111)/2

   M122_1 <- sum(Ubb*Ybb*YT)/2

   M122_2 <- sum(UVbb*Ybb)

   M122_3 <- (sum(Ubb*Two(Vbb))-M112)/2

   M123_1 <- sum(Ubb*Two(Ybb)*YT)

   M123_2 <- sum(UVbb*Two(YT))

   M123_3 <- sum(UVbb*Ybb*YT)/2

   M123_4 <- sum(Ubb*Two(Vbb)*Ybb)

   M123_5 <- (sum(Ubb*Three(Vbb))-M113)/2

   M221_1 <- sum(Xbb*XT*Vbb)/2

   M221_2 <- sum(UVbb*Xbb)

   M221_3 <- (sum(Two(Ubb)*Vbb)-M211)/2

   M222_1 <- sum(Ubb*Xbb*Ybb*YT)

   M222_2 <- sum(Two(Ubb)*Ybb*YT)/2

   M222_3 <- sum(Xbb*XT*Ybb*Vbb)

   M222_4 <- sum(UVbb*Xbb*Ybb)

   M222_5 <- sum(Two(Ubb)*Vbb*YT)

   M222_6 <- sum(Xbb*XT*Two(Vbb))/2

   M222_7 <- sum(Ubb*Xbb*Two(Vbb))

   M222_8 <- (sum(Two(Ubb)*Two(Vbb))-M212)/2

   M222_9 <- sum(UVbb*Xbb*YT)

   M321_1 <- sum(Two(Xbb)*XT*Vbb)

   M321_2 <- sum(UVbb*Two(XT))

   M321_3 <- sum(UVbb*Xbb*XT)/2

   M321_4 <- sum(Two(Ubb)*XT*Vbb)

   M321_5 <- (sum(Three(Ubb)*Vbb)-M311)/2

   M131 <- sum(Three(PQ))

   M141 <- sum(Four(PQ))

   M1321<-M1322<-M1323<-M1324<-M1325<-M2311<-M2312<-M2313<-M2314<-M2315<-0
   for(i in 1:L1){
      for(j in 1:L1){
         if(i!=j){
            P_F <- colSums(P[,i]*P[,j]*P[,-c(i,j)])
            P_F1 <- colSums(P[,i]*P[,j]*(1-P[,-c(i,j)]))
            P_F2 <- colSums((1-P[,i])*(1-P[,j])*P[,-c(i,j)])
            P_F3 <- colSums((1-P[,i])*P[,j]*P[,-c(i,j)])
            Q_F <- colSums(Q[i,]*Q[j,]*t(Q[-c(i,j),]))
            Q_F1 <- colSums(Q[i,]*Q[j,]*t(1-Q[-c(i,j),]))
            Q_F2 <- colSums((1-Q[i,])*(1-Q[j,])*t(Q[-c(i,j),]))
            Q_F3 <- colSums((1-Q[i,])*Q[j,]*t(Q[-c(i,j),]))

            M1321<-M1321+sum(P_F*Q_F1*Q_F2)
            M1322<-M1322+sum(P_F*Q_F*Q_F2)
            M1323<-M1323+sum(P_F*Q_F1*Q_F3)
            M1324<-M1324+sum(P_F*Q_F*Q_F1)
            M1325<-M1325+sum(P_F*Q_F*(Q_F-1)/2)
            M2311<-M2311+sum(P_F1*P_F2*Q_F)
            M2312<-M2312+sum(P_F*P_F2*Q_F)
            M2313<-M2313+sum(P_F1*P_F3*Q_F)
            M2314<-M2314+sum(P_F*P_F3*Q_F)
            M2315<-M2315+sum((P_F*(P_F-1)/2)*Q_F)
         }
      }
   }
   M1321<-M1321/2
   M1322<-M1322/2
   M1323<-M1323/2
   M1324<-M1324/2
   M1325<-M1325/6
   M2311<-M2311/2
   M2312<-M2312/2
   M2313<-M2313/2
   M2314<-M2314/2
   M2315<-M2315/6
   motif<-c(M111, M112, M113, M114, M211, M212, M213, M311, M312, M411,
            M121, M122_1, M122_2, M122_3, M123_1, M123_2, M123_3,
            M123_4, M123_5, M221_1, M221_2, M221_3, M222_1, M222_2,
            M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_9, M321_1,
            M321_2, M321_3, M321_4, M321_5, M131, M1321, M1322,
            M1323, M1324, M1325, M2311, M2312, M2313, M2314, M2315, M141)

   motif_names<-c("M111", "M112","M113", "M114", "M211", "M212", "M213", "M311", "M312", "M411",
                  "M121", "M122-1", "M122-2", "M122-3", "M123-1", "M123-2", "M123-3",
                  "M123-4", "M123-5", "M221-1", "M221-2", "M221-3", "M222-1", "M222-2",
                  "M222-3", "M222-4", "M222-5", "M222-6", "M222-7", "M222-8", "M222-9", "M321-1",
                  "M321-2", "M321-3", "M321-4", "M321-5", "M131", "M132-1", "M132-2",
                  "M132-3", "M132-4", "M132-5", "M231-1", "M231-2", "M231-3", "M231-4", "M231-5", "M141")
    if(weighted){

      ObW <- colSums(PW)
      RbW <- rowSums(QW)

      B_A2 <- apply(PW,2,function(x){ (sum(x!=0)-1)*sum(x) })
      B_A3 <- apply(PW,2,function(x){ n <- sum(x!=0); (n-1)*(n-2)*sum(x)/2 })
      B_A4 <- apply(PW,2,function(x){ n <- sum(x!=0); (n-1)*(n-2)*(n-3)*sum(x)/6 })

      B_C2 <- apply(QW,1,function(x){ (sum(x!=0)-1)*sum(x) })
      B_C3 <- apply(QW,1,function(x){ n <- sum(x!=0); (n-1)*(n-2)*sum(x)/2 })
      B_C4 <- apply(QW,1,function(x){ n <- sum(x!=0); (n-1)*(n-2)*(n-3)*sum(x)/6 })


      M111 <- sum(ObW * Rb + RbW * Ob)

      M112 <- sum(ObW * Two(Rb) + B_C2 * Ob)

      M113 <- sum(ObW * Three(Rb) + B_C3 * Ob)

      M114 <- sum(ObW * Four(Rb) + B_C4 * Ob)

      M211 <- sum(B_A2 * Rb + RbW * Two(Ob))

      M212 <- sum(B_A2 * Two(Rb) + B_C2 * Two(Ob))

      M213 <- sum(B_A2 * Three(Rb) + B_C3 * Two(Ob))

      M311 <- sum(B_A3 * Rb + RbW * Three(Ob))

      M312 <- sum(B_A3 * Two(Rb) + B_C2 * Three(Ob))

      M411 <- sum(B_A4 * Rb + RbW * Four(Ob))

      All_add <- function(a,b){ return((a!=0)*(b!=0)*(a+b)) }

      A_2BW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum(All_add(x,y)) })}))
      C_2BW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum(All_add(x,y)) })}))

      B_AW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum((x!=0)*(y==0)*(x+y)) })}))
      B_ATW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum((x==0)*(y!=0)*(x+y)) })}))
      B_CW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum((x!=0)*(y==0)*(x+y)) })}))
      B_CTW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum((x==0)*(y!=0)*(x+y)) })}))

      M121 <- (sum(A_2BW*Vbb + C_2BW*Ubb) - 2*M111)/2

      M122_1 <- sum(A_2BW*Ybb*YT + (B_CW *YT + B_CTW *Ybb) *Ubb)/2

      M122_2 <- sum(A_2BW*Vbb*YT + (C_2BW *YT + B_CTW *Vbb) *Ubb)

      M122_3 <- (sum(A_2BW*Vbb*(Vbb-1)/2 + C_2BW *(Vbb-1) *Ubb) - 2*M112)/2

      M123_1 <- sum(A_2BW*Ybb*(Ybb-1)*YT/2 + (B_CW *(Ybb-1)*YT + B_CTW *Ybb*(Ybb-1)/2) *Ubb)

      M123_2 <- sum(A_2BW*Vbb*YT*(YT-1)/2 + (C_2BW *YT*(YT-1)/2 + B_CTW *(YT-1)*Vbb) *Ubb)

      M123_3 <- sum(A_2BW*Vbb*Ybb*YT + (C_2BW *Ybb*YT + B_CW * Vbb * YT + B_CTW *Vbb * Ybb) *Ubb)/2

      M123_4 <- sum(A_2BW*Vbb*(Vbb-1)*Ybb/2 + (C_2BW *(Vbb-1)*Ybb + B_CW * Vbb * (Vbb-1)/2) *Ubb)

      M123_5 <- (sum(A_2BW*Vbb*(Vbb-1)*(Vbb-2)/6 + (C_2BW *(Vbb-1)*(Vbb-2)/2 ) *Ubb) -2 * M113)/2

      M221_1 <- sum((B_AW *XT + B_ATW *Xbb) *Vbb + C_2BW*Xbb*XT  )/2

      M221_2 <- sum((B_AW *Ubb + A_2BW *Xbb) *Vbb + C_2BW*Xbb*Ubb )

      M221_3 <- (sum(A_2BW*(Ubb-1)*Vbb + C_2BW *Ubb*(Ubb-1)/2) - 2*M211)/2

      M222_1 <- sum((B_AW * Ubb + A_2BW * Xbb) *Ybb*YT +(B_CW*YT+B_CTW*Ybb) *Xbb*Ubb)

      M222_2 <- sum(A_2BW *(Ubb-1) *Ybb*YT +(B_CW*YT+B_CTW*Ybb) *Ubb*(Ubb-1)/2)/2

      M222_3 <- sum((B_AW * XT + B_ATW * Xbb) *Ybb*Vbb +(B_CW*Vbb + C_2BW*Ybb) *Xbb*XT)

      M222_4 <- sum((B_AW * Ubb + A_2BW * Xbb) *Ybb*Vbb +(B_CW*Vbb + C_2BW*Ybb) *Xbb*Ubb)

      M222_5 <- sum((A_2BW * (Ubb-1)*Vbb*YT +(C_2BW*YT + B_CTW* Vbb) *Ubb*(Ubb-1)/2))

      M222_6 <- sum((B_AW * XT + B_ATW * Xbb) * Vbb*(Vbb-1)/2 + C_2BW* (Vbb-1) *Xbb*XT)/2

      M222_7 <- sum((B_AW * Ubb + A_2BW * Xbb) *Vbb*(Vbb-1)/2 + C_2BW*(Vbb-1) *Xbb*Ubb)

      M222_8 <- (sum(A_2BW *(Ubb-1) *Vbb*(Vbb-1)/2 + C_2BW*(Vbb-1) *Ubb*(Ubb-1)/2) -2*M212)/2

      M222_9 <- sum((B_AW * Ubb + A_2BW * Xbb) *Vbb*YT +(C_2BW*YT + B_CTW*Vbb) *Xbb*Ubb)

      M321_1 <- sum((B_AW *(Xbb-1)*XT + B_ATW *Xbb*(Xbb-1)/2) *Vbb + C_2BW *Xbb*(Xbb-1)*XT/2)

      M321_2 <- sum((A_2BW *XT*(XT-1)/2 + B_ATW *(XT-1)*Ubb) *Vbb + C_2BW *Ubb*XT*(XT-1)/2)

      M321_3 <- sum((B_AW *Ubb*XT + A_2BW * Xbb*XT + B_ATW *Xbb*Ubb) *Vbb + C_2BW *Xbb*Ubb*XT)/2

      M321_4 <- sum((A_2BW *(Ubb-1)*XT + B_ATW *Ubb*(Ubb-1)/2) *Vbb + C_2BW *Ubb*(Ubb-1)*XT/2)

      M321_5 <- (sum((A_2BW *(Ubb-1)*(Ubb-2)/2) *Vbb + C_2BW *Ubb*(Ubb-1)*(Ubb-2)/6) -2*M311)/2


      PQW <- t(apply(PW, 1, function(x){apply(QW,2, function(y){sum(All_add(x,y))})}))

      M131 <- sum(PQW*(PQ)*(PQ-2)/2)

      M141 <- sum(PQW*(PQ-1)*(PQ-2)*(PQ-3)/6)

      M1321<-M1322<-M1323<-M1324<-M1325<-M2311<-M2312<-M2313<-M2314<-M2315<-0

      P0 <- (PW==0)*1; Q0 <- (Q==0)*1

      for(i in 1:L1){
         for(j in 1:L1){
            if(i!=j){
               P_F <- P[,i]*P[,j]*P[,-c(i,j)]
               P_F1 <- P[,i]*P[,j]*P0[,-c(i,j)]
               P_F2<- P0[,i]*P0[,j]*P[,-c(i,j)]
               P_F3 <- P0[,i]*P[,j]*P[,-c(i,j)]
               Q_F <- Q[i,]*Q[j,]*t(Q[-c(i,j),])
               Q_F1 <- Q[i,]*Q[j,]*t(Q0[-c(i,j),])
               Q_F2 <- Q0[i,]*Q0[j,]*t(Q[-c(i,j),])
               Q_F3 <- Q0[i,]*Q[j,]*t(Q[-c(i,j),])

               P_FCS <- colSums(P_F)
               P_F1CS <- colSums(P_F1)
               P_F2CS<- colSums(P_F2)
               P_F3CS <- colSums(P_F3)
               Q_FCS <- colSums(Q_F)
               Q_F1CS <- colSums(Q_F1)
               Q_F2CS <- colSums(Q_F2)
               Q_F3CS <- colSums(Q_F3)

               PW_ijk <- PW[,i]+PW[,j]+PW[,-c(i,j)]
               QW_ijk <- QW[i,]+QW[j,]+t(QW[-c(i,j),])

               PW_F <- colSums(PW_ijk*P_F)
               PW_F1 <- colSums(PW_ijk*P_F1)
               PW_F2 <- colSums(PW_ijk*P_F2)
               PW_F3 <- colSums(PW_ijk*P_F3)
               QW_F <- colSums(QW_ijk*Q_F)
               QW_F1 <- colSums(QW_ijk*Q_F1)
               QW_F2 <- colSums(QW_ijk*Q_F2)
               QW_F3 <- colSums(QW_ijk*Q_F3)

               M1321<-M1321+sum(PW_F*Q_F1CS*Q_F2CS + (QW_F1*Q_F2CS +QW_F2*Q_F1CS)*P_FCS)
               M1322<-M1322+sum(PW_F*Q_FCS*Q_F2CS + (QW_F*Q_F2CS +QW_F2*Q_FCS)*P_FCS)
               M1323<-M1323+sum(PW_F*Q_F1CS*Q_F3CS + (QW_F1*Q_F3CS +QW_F3*Q_F1CS)*P_FCS)
               M1324<-M1324+sum(PW_F*Q_FCS*Q_F1CS + (QW_F*Q_F1CS +QW_F1*Q_FCS)*P_FCS)
               M1325<-M1325+sum(PW_F*Q_FCS*(Q_FCS-1)/2 + QW_F*(Q_FCS-1)*P_FCS)
               M2311<-M2311+sum((PW_F1*P_F2CS+PW_F2*P_F1CS)*Q_FCS + QW_F*P_F1CS*P_F2CS)
               M2312<-M2312+sum((PW_F*P_F2CS+PW_F2*P_FCS)*Q_FCS + QW_F*P_FCS*P_F2CS)
               M2313<-M2313+sum((PW_F1*P_F3CS+PW_F3*P_F1CS)*Q_FCS + QW_F*P_F1CS*P_F3CS)
               M2314<-M2314+sum((PW_F*P_F3CS+PW_F3*P_FCS)*Q_FCS + QW_F*P_FCS*P_F3CS)
               M2315<-M2315+sum(PW_F*(P_FCS-1)*Q_FCS + QW_F*P_FCS*(P_FCS-1)/2)
            }
         }
      }
      M1321<-M1321/2
      M1322<-M1322/2
      M1323<-M1323/2
      M1324<-M1324/2
      M1325<-M1325/6
      M2311<-M2311/2
      M2312<-M2312/2
      M2313<-M2313/2
      M2314<-M2314/2
      M2315<-M2315/6

      edge_php <- c(2,3,4,5,3,4,5,4,5,5,4,4,5,6,5,6,6,7,8,4,5,6,5,6,5,6,7,6,7,8,6,5,6,6,7,8,6,6,7,7,8,9,6,7,7,8,9,8)

      motif_weighted<-c(M111, M112, M113, M114, M211, M212, M213, M311, M312, M411,
                        M121, M122_1, M122_2, M122_3, M123_1, M123_2, M123_3,
                        M123_4, M123_5, M221_1, M221_2, M221_3, M222_1, M222_2,
                        M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_9, M321_1,
                        M321_2, M321_3, M321_4, M321_5, M131, M1321, M1322,
                        M1323, M1324, M1325, M2311, M2312, M2313, M2314, M2315, M141)
      mean_weighted <- motif_weighted/(motif*edge_php)

      mean_weighted <- replace(mean_weighted,which(is.nan(mean_weighted)),NA)

      return(data.frame(motif_name=motif_names,count=motif,mean_weight=mean_weighted ))
   }
   else{
      return(data.frame(motif_name=motif_names,count=motif))
   }

}
