#' Count the roles of connector nodes defined by interconnection motifs
#'
#' Count the roles of connector nodes defined by interconnection motifs in a tripartite network.
#'
#' @param network.or.subnet_mat1 An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.
#' @param subnet_mat2 A matrix representing one subnetwork.
#' @param weighted Logical. Default to FALSE. If TRUE, a weighted measure is provided. See details.

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
#' <br>For weighted tripartite networks, the mean weight of the motif occurrence (i.e., a motif occurrence isomorphic to a particular motif form) is provided for a given node with a given role, following Mora et al. (2018) and Simmons et al. (2019).
#'
#' @import igraph
#' @export
#'
#' @return
#' For binary networks, return a matrix with elements representing the number of times each connector node plays for each unique role within interconnection motifs; for weighted networks, the matrix element represents the mean weight of the motif occurrences where the node exists.
#'
#'
#' @srrstats {G1.1} The algorithm is the first implementation of a novel algorithm.
#'
#' @references
#' #' Mora, B.B., Cirtwill, A.R. and Stouffer, D.B. (2018). pymfinder: a tool for the motif analysis of binary and quantitative complex networks. bioRxiv, 364703.
#'
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_toy_net(11,15,16,0.2)
#' icmotif_role(Net)
#'
#' ## empirical network
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' icmotif_role(Net)
#' set.seed(13)
#' library(igraph)
#' E(Net)$weight<-runif(length(E(Net)),0.1,1)#random weights assigned
#' icmotif_role(Net, weighted=T)
#'
#'
#' ##input as binary matrices, with row names.
#' set.seed(12)
#' md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
#' md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
#' icmotif_role(md1,md2)
#'
#'
#'##input as weighted matrices,with row numbers as row names.
#' set.seed(12)
#' mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
#' mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' icmotif_role(mdw1,mdw2,weighted=T)
#'

icmotif_role<-function(network.or.subnet_mat1, subnet_mat2=NULL, weighted =F){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network <- adjust_net(network.or.subnet_mat1)
      PQ <- as.matrix(network[])
      P <- PQ[(V(network)$level)==0,(V(network)$level)==1]
      Q <- PQ[(V(network)$level)==1,(V(network)$level)==2]
      spe <- V(network.or.subnet_mat1)$name[V(network.or.subnet_mat1)$level==1]
      role <- matrix(0,length(spe),70)
      rownames(role) <- spe
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
      mat1 <- network.or.subnet_mat1
      mat2 <- subnet_mat2

      if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
         rownames(mat1) <- paste0("mid_spe",seq=1:nrow(mat1))
         rownames(mat2) <- paste0("mid_spe",seq=1:nrow(mat2))
         matrow <- unique(c(rownames(mat1),rownames(mat2)))
      }
      if(nrow(mat1)!=nrow(mat2))
         message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
      if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
         matrow <- unique(c(rownames(mat1),rownames(mat2)))
      else
         stop("Make sure matrices either have no row names or have full row names. No NA!!!")
      mat_1 <- matrix(0,length(matrow),ncol(mat1))
      rownames(mat_1) <- matrow
      mat_1[rownames(mat1),] <- mat1

      mat_2<-matrix(0,length(matrow),ncol(mat2))
      rownames(mat_2) <- matrow
      mat_2[rownames(mat2),] <- mat2

      P <- t(mat_1)
      Q <- mat_2

      logi <- (apply(P,2,sum)*apply(Q,1,sum))!=0
      P <- P[,logi]
      Q <- Q[logi,]
      spe <- rownames(mat_1)[rownames(mat_1)==rownames(mat_2)]
      role <- matrix(0,length(spe),70)
      rownames(role) <- spe
   }
   else
      stop("Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   L1 <- ncol(P)
   if(L1 < 4L)
      stop("Error: please input a large 'number of interconnecting species >=4' network data!!!")
   PW <- P; QW <- Q
   P <- (P!=0)*1; Q <- (Q!=0)*1


   Ob <- colSums(P)
   Rb <- rowSums(Q)
   Ubb <- t(P) %*% P
   Xbb <- t(P) %*% (1 - P)
   XT <- t(Xbb)
   Vbb <- Q %*% t(Q)
   Ybb <- Q %*% (1 - t(Q))
   YT <- t(Ybb)
   UVbb <- Ubb * Vbb
   Two <-function(a){ return(a*(a-1)/2) }
   Three <- function(a){ return(a*(a-1)*(a-2)/6)}
   Four <- function(a){ return(a*(a-1)*(a-2)*(a-3)/24)}
   BB_role <- NULL

   M111 <- Ob*Rb
   BB_role <- cbind(BB_role,M111)

   M112 <- Ob*Two(Rb)
   BB_role <- cbind(BB_role,M112)

   M113 <- Ob*Three(Rb)
   BB_role <- cbind(BB_role,M113)

   M114<-Ob*Four(Rb)
   BB_role <- cbind(BB_role,M114)

   M211 <- Two(Ob)*Rb
   BB_role <- cbind(BB_role,M211)

   M212 <- Two(Ob)*Two(Rb)
   BB_role <- cbind(BB_role,M212)

   M213 <- Two(Ob)*Three(Rb)
   BB_role <- cbind(BB_role,M213)

   M311 <- Three(Ob)*Rb
   BB_role <- cbind(BB_role,M311)

   M312 <- Three(Ob)*Two(Rb)
   BB_role <- cbind(BB_role,M312)

   M411 <- Four(Ob)*Rb
   BB_role <- cbind(BB_role,M411)

   # M121#
   BB_role <- cbind(BB_role,rowSums(UVbb)-M111)

   M122_1 <- Ubb*Ybb*YT
   BB_role <- cbind(BB_role,rowSums(M122_1))

   M122_2 <- UVbb*Ybb
   BB_role <- cbind(BB_role,colSums(M122_2),rowSums(M122_2))

   M122_3 <- Ubb*Two(Vbb)
   BB_role <- cbind(BB_role,rowSums(M122_3)-M112)

   M123_1 <- Ubb*Two(Ybb)*YT
   BB_role <- cbind(BB_role,rowSums(M123_1),colSums(M123_1))

   M123_2 <- UVbb*Two(YT)
   BB_role <- cbind(BB_role,rowSums(M123_2),colSums(M123_2))

   M123_3 <- UVbb*Ybb*YT
   BB_role <- cbind(BB_role,rowSums(M123_3))

   M123_4 <- Ubb*Two(Vbb)*Ybb
   BB_role <- cbind(BB_role,rowSums(M123_4),colSums(M123_4))

   M123_5 <- Ubb*Three(Vbb)
   BB_role <- cbind(BB_role,rowSums(M123_5)-M113)

   M221_1 <- Xbb*XT*Vbb
   BB_role <- cbind(BB_role,rowSums(M221_1))

   M221_2 <- UVbb*Xbb
   BB_role <- cbind(BB_role,rowSums(M221_2),colSums(M221_2))

   M221_3 <- Two(Ubb)*Vbb
   BB_role <- cbind(BB_role,rowSums(M221_3)-M211)

   M222_1 <- Ubb*Xbb*Ybb*YT
   BB_role <- cbind(BB_role,rowSums(M222_1),colSums(M222_1))

   M222_2 <- Two(Ubb)*Ybb*YT
   BB_role <- cbind(BB_role,rowSums(M222_2))

   M222_3 <- Xbb*XT*Ybb*Vbb
   BB_role <- cbind(BB_role,rowSums(M222_3),colSums(M222_3))

   M222_4 <- UVbb*Xbb*Ybb
   BB_role <- cbind(BB_role,rowSums(M222_4),colSums(M222_4))

   M222_5 <- Two(Ubb)*Vbb*YT
   BB_role <- cbind(BB_role,rowSums(M222_5),colSums(M222_5))

   M222_6 <- Xbb*XT*Two(Vbb)
   BB_role <- cbind(BB_role,rowSums(M222_6))

   M222_7 <- Ubb*Xbb*Two(Vbb)
   BB_role <- cbind(BB_role,rowSums(M222_7),colSums(M222_7))

   M222_8 <- Two(Ubb)*Two(Vbb)
   BB_role <- cbind(BB_role,rowSums(M222_8)-M212)

   M222_9 <- UVbb*Xbb*YT
   BB_role <- cbind(BB_role,rowSums(M222_9),colSums(M222_9))

   M321_1 <- Two(Xbb)*XT*Vbb
   BB_role <- cbind(BB_role,rowSums(M321_1),colSums(M321_1))

   M321_2 <- UVbb*Two(XT)
   BB_role <- cbind(BB_role,rowSums(M321_2),colSums(M321_2))

   M321_3 <- UVbb*Xbb*XT
   BB_role <- cbind(BB_role,rowSums(M321_3))

   M321_4 <- Two(Ubb)*XT*Vbb
   BB_role <- cbind(BB_role,rowSums(M321_4),colSums(M321_4))

   M321_5 <- Three(Ubb)*Vbb
   BB_role <- cbind(BB_role,rowSums(M321_5)-M311)


   ####################
   M131<-M1321<-M1321_1<-M1322<-M1322_1<-M1323<-M1323_1<-M1324<-rep(0,L1)
   M1324_1<-M1325<-M2311<-M2311_1<-M2312<-M2312_1<-M2313<-M2313_1<-rep(0,L1)
   M2314=M2314_1=M2315=M141<-rep(0,L1)
   if(L1>4){
      for(i in 1:L1){
         for(j in 1:L1){
            if(i!=j){
               P_F<-colSums(P[,i]*P[,j]*(P[,-c(i,j)]))
               P_F1<-colSums(P[,i]*P[,j]*(1-P[,-c(i,j)]))
               P_F2<-colSums((1-P[,i])*(1-P[,j])*(P[,-c(i,j)]))
               P_F3<-colSums((1-P[,i])*(P[,j])*(P[,-c(i,j)]))
               Q_F<-colSums(Q[i,]*Q[j,]*t(Q[-c(i,j),]))
               Q_F1<-colSums(Q[i,]*Q[j,]*t(1-Q[-c(i,j),]))
               Q_F2<-colSums((1-Q[i,])*(1-Q[j,])*t(Q[-c(i,j),]))
               Q_F3<-colSums((1-Q[i,])*(Q[j,])*t(Q[-c(i,j),]))
               ########
               role_value<-P_F*Q_F
               if(sum(role_value)!=0){
                  M131[c(i,j)]<-M131[c(i,j)]+sum(role_value)
                  M131[-c(i,j)]<-M131[-c(i,j)]+role_value
               }
               role_value1<-P_F*Q_F1*Q_F2
               if(sum(role_value1)!=0){
                  M1321[c(i,j)]<-M1321[c(i,j)]+sum(role_value1)
                  M1321_1[-c(i,j)]<-M1321_1[-c(i,j)]+role_value1
               }
               role_value2<-P_F*Q_F*Q_F2
               if(sum(role_value2)!=0){
                  M1322[c(i,j)]<-M1322[c(i,j)]+sum(role_value2)
                  M1322_1[-c(i,j)]<-M1322_1[-c(i,j)]+role_value2
               }
               role_value3<-P_F*Q_F1*Q_F3
               if(sum(role_value3)!=0){
                  M1323[i]<-M1323[i]+sum(role_value3)
                  M1323[-c(i,j)]<-M1323[-c(i,j)]+role_value3
                  M1323_1[j]<-M1323_1[j]+sum(role_value3)
               }
               role_value4<-P_F*Q_F*Q_F1
               if(sum(role_value4)!=0){
                  M1324[c(i,j)]<-M1324[c(i,j)]+sum(role_value4)
                  M1324_1[-c(i,j)]<-M1324_1[-c(i,j)]+role_value4
               }
               role_value5<-P_F*Q_F*(Q_F-1)/2
               if(sum(role_value5)!=0){
                  M1325[c(i,j)]<-M1325[c(i,j)]+sum(role_value5)
                  M1325[-c(i,j)]<-M1325[-c(i,j)]+role_value5
               }
               role_value6<-P_F1*P_F2*Q_F
               if(sum(role_value6)!=0){
                  M2311[c(i,j)]<-M2311[c(i,j)]+sum(role_value6)
                  M2311_1[-c(i,j)]<-M2311_1[-c(i,j)]+role_value6
               }
               role_value7<-P_F*P_F2*Q_F
               if(sum(role_value7)!=0){
                  M2312_1[c(i,j)]<-M2312_1[c(i,j)]+sum(role_value7)
                  M2312[-c(i,j)]<-M2312[-c(i,j)]+(role_value7)
               }
               role_value8<-P_F1*P_F3*Q_F
               if(sum(role_value8)!=0){
                  M2313[i]<-M2313[i]+sum(role_value8)
                  M2313[-c(i,j)]<-M2313[-c(i,j)]+role_value8
                  M2313_1[j]<-M2313_1[j]+sum(role_value8)
               }
               role_value9<-P_F*P_F3*Q_F
               if(sum(role_value9)!=0){
                  M2314[i]<-M2314[i]+sum(role_value9)
                  M2314_1[j]<-M2314_1[j]+sum(role_value9)
                  M2314_1[-c(i,j)]<-M2314_1[-c(i,j)]+role_value9
               }
               role_value10<-P_F*(P_F-1)*Q_F/2
               if(sum(role_value10)!=0){
                  M2315[c(i,j)]<-M2315[c(i,j)]+sum(role_value10)
                  M2315[-c(i,j)]<-M2315[-c(i,j)]+role_value10
               }
            }
            if(i<j){
               for(k in (j):L1){
                  if(j<k){
                     role_value11<-((P[,i]*P[,j]*P[,k])%*%P[,-c(i,j,k)])*((Q[i,]*Q[j,]*Q[k,])%*%t(Q[-c(i,j,k),]))
                     if(sum(role_value11)!=0){
                        M141[c(i,j,k)]<-M141[c(i,j,k)]+sum(role_value11)
                        M141[-c(i,j,k)]<-M141[-c(i,j,k)]+role_value11
                     }
                  }
               }
            }
         }
      }
      M141<-M141/4
   }
   else{
      for(i in 1:4){
         for(j in 1:4){
            if(i!=j){
               P_F<-colSums(P[,i]*P[,j]*(P[,-c(i,j)]))
               P_F1<-colSums(P[,i]*P[,j]*(1-P[,-c(i,j)]))
               P_F2<-colSums((1-P[,i])*(1-P[,j])*(P[,-c(i,j)]))
               P_F3<-colSums((1-P[,i])*(P[,j])*(P[,-c(i,j)]))
               Q_F<-colSums(Q[i,]*Q[j,]*t(Q[-c(i,j),]))
               Q_F1<-colSums(Q[i,]*Q[j,]*t(1-Q[-c(i,j),]))
               Q_F2<-colSums((1-Q[i,])*(1-Q[j,])*t(Q[-c(i,j),]))
               Q_F3<-colSums((1-Q[i,])*(Q[j,])*t(Q[-c(i,j),]))
               ########
               role_value<-P_F*Q_F
               if(sum(role_value)!=0){
                  M131[c(i,j)]<-M131[c(i,j)]+sum(role_value)
                  M131[-c(i,j)]<-M131[-c(i,j)]+role_value
               }
               role_value1<-P_F*Q_F1*Q_F2
               if(sum(role_value1)!=0){
                  M1321[c(i,j)]<-M1321[c(i,j)]+sum(role_value1)
                  M1321_1[-c(i,j)]<-M1321_1[-c(i,j)]+role_value1
               }
               role_value2<-P_F*Q_F*Q_F2
               if(sum(role_value2)!=0){
                  M1322[c(i,j)]<-M1322[c(i,j)]+sum(role_value2)
                  M1322_1[-c(i,j)]<-M1322_1[-c(i,j)]+role_value2
               }
               role_value3<-P_F*Q_F1*Q_F3
               if(sum(role_value3)!=0){
                  M1323[i]<-M1323[i]+sum(role_value3)
                  M1323[-c(i,j)]<-M1323[-c(i,j)]+role_value3
                  M1323_1[j]<-M1323_1[j]+sum(role_value3)
               }
               role_value4<-P_F*Q_F*Q_F1
               if(sum(role_value4)!=0){
                  M1324[c(i,j)]<-M1324[c(i,j)]+sum(role_value4)
                  M1324_1[-c(i,j)]<-M1324_1[-c(i,j)]+role_value4
               }
               role_value5<-P_F*Q_F*(Q_F-1)/2
               if(sum(role_value5)!=0){
                  M1325[c(i,j)]<-M1325[c(i,j)]+sum(role_value5)
                  M1325[-c(i,j)]<-M1325[-c(i,j)]+role_value5
               }
               role_value6<-P_F1*P_F2*Q_F
               if(sum(role_value6)!=0){
                  M2311[c(i,j)]<-M2311[c(i,j)]+sum(role_value6)
                  M2311_1[-c(i,j)]<-M2311_1[-c(i,j)]+role_value6
               }
               role_value7<-P_F*P_F2*Q_F
               if(sum(role_value7)!=0){
                  M2312_1[c(i,j)]<-M2312_1[c(i,j)]+sum(role_value7)
                  M2312[-c(i,j)]<-M2312[-c(i,j)]+(role_value7)
               }
               role_value8<-P_F1*P_F3*Q_F
               if(sum(role_value8)!=0){
                  M2313[i]<-M2313[i]+sum(role_value8)
                  M2313[-c(i,j)]<-M2313[-c(i,j)]+role_value8
                  M2313_1[j]<-M2313_1[j]+sum(role_value8)
               }
               role_value9<-P_F*P_F3*Q_F
               if(sum(role_value9)!=0){
                  M2314[i]<-M2314[i]+sum(role_value9)
                  M2314_1[j]<-M2314_1[j]+sum(role_value9)
                  M2314_1[-c(i,j)]<-M2314_1[-c(i,j)]+role_value9
               }
               role_value10<-P_F*(P_F-1)*Q_F/2
               if(sum(role_value10)!=0){
                  M2315[c(i,j)]<-M2315[c(i,j)]+sum(role_value10)
                  M2315[-c(i,j)]<-M2315[-c(i,j)]+role_value10
               }
            }
            M141<-sum(P[,1]*P[,2]*P[,3]*P[,4])*sum(Q[1,]*Q[2,]*Q[3,]*Q[4,])
         }
      }
   }
   BB_role<-cbind(BB_role,M131/6,M1321/2,M1321_1/2,M1322_1/2,M1322/2,M1323/2,M1323_1/2,M1324/2,M1324_1/2,M1325/6,M2311/2,M2311_1/2,M2312/2,M2312_1/2,M2313/2,M2313_1/2,M2314/2,M2314_1/2,M2315/6,M141)
   colnames(BB_role)<-paste0("role",c(1:70))
   role[rownames(BB_role),]<-BB_role
   colnames(role)<-paste0("role",c(1:70))


   if(weighted){
      ObW <- colSums(PW)
      RbW <- rowSums(QW)

      B_A2 <- apply(PW,2,function(x){ (sum(x!=0)-1)*sum(x) })
      B_A3 <- apply(PW,2,function(x){ n <- sum(x!=0); (n-1)*(n-2)*sum(x)/2 })
      B_A4 <- apply(PW,2,function(x){ n <- sum(x!=0); (n-1)*(n-2)*(n-3)*sum(x)/6 })

      B_C2 <- apply(QW,1,function(x){ (sum(x!=0)-1)*sum(x) })
      B_C3 <- apply(QW,1,function(x){ n <- sum(x!=0); (n-1)*(n-2)*sum(x)/2 })
      B_C4 <- apply(QW,1,function(x){ n <- sum(x!=0); (n-1)*(n-2)*(n-3)*sum(x)/6 })

      BB_role <- NULL

      M111 <- (ObW * Rb + RbW * Ob)/2
      BB_role <- cbind(BB_role,M111)

      M112 <- (ObW * Two(Rb) + B_C2 * Ob)/3
      BB_role <- cbind(BB_role,M112)


      M113 <- (ObW * Three(Rb) + B_C3 * Ob)/4
      BB_role <- cbind(BB_role,M113)

      M114 <- (ObW * Four(Rb) + B_C4 * Ob)/5
      BB_role <- cbind(BB_role,M114)

      M211 <- (B_A2 * Rb + RbW * Two(Ob))/3
      BB_role <- cbind(BB_role,M211)

      M212 <- (B_A2 * Two(Rb) + B_C2 * Two(Ob))/4
      BB_role <- cbind(BB_role,M212)

      M213 <- (B_A2 * Three(Rb) + B_C3 * Two(Ob))/5
      BB_role <- cbind(BB_role,M213)

      M311 <- (B_A3 * Rb + RbW * Three(Ob))/4
      BB_role <- cbind(BB_role,M311)

      M312 <- (B_A3 * Two(Rb) + B_C2 * Three(Ob))/5
      BB_role <- cbind(BB_role,M312)

      M411 <- (B_A4 * Rb + RbW * Four(Ob))/5
      BB_role <- cbind(BB_role,M411)

      All_add <- function(a,b){ return((a!=0)*(b!=0)*(a+b)) }

      A_2BW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum(All_add(x,y)) })}))
      C_2BW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum(All_add(x,y)) })}))

      B_AW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum((x!=0)*(y==0)*(x+y)) })}))
      B_ATW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum((x==0)*(y!=0)*(x+y)) })}))
      B_CW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum((x!=0)*(y==0)*(x+y)) })}))
      B_CTW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum((x==0)*(y!=0)*(x+y)) })}))

      M121 <- rowSums(A_2BW*Vbb + C_2BW*Ubb)/4 - M111
      BB_role <- cbind(BB_role,M121)

      M122_1 <- rowSums(A_2BW*Ybb*YT + (B_CW *YT + B_CTW *Ybb) *Ubb)/4
      BB_role <- cbind(BB_role,M122_1)

      M122_2 <- (A_2BW*Vbb*YT + (C_2BW *YT + B_CTW *Vbb) *Ubb)/5
      BB_role <- cbind(BB_role,rowSums(M122_2),colSums(M122_2))

      M122_3 <- rowSums(A_2BW*Vbb*(Vbb-1)/2 + C_2BW *(Vbb-1) *Ubb)/6 - M112
      BB_role <- cbind(BB_role,M122_3)


      M123_1 <- (A_2BW*Ybb*(Ybb-1)*YT/2 + (B_CW *(Ybb-1)*YT + B_CTW *Ybb*(Ybb-1)/2) *Ubb)/5
      BB_role <- cbind(BB_role,rowSums(M123_1),colSums(M123_1))

      M123_2 <- (A_2BW*Vbb*YT*(YT-1)/2 + (C_2BW *YT*(YT-1)/2 + B_CTW *(YT-1)*Vbb) *Ubb)/6
      BB_role <- cbind(BB_role,rowSums(M123_2),colSums(M123_2))

      M123_3 <- rowSums(A_2BW*Vbb*Ybb*YT + (C_2BW *Ybb*YT + B_CW * Vbb * YT + B_CTW *Vbb * Ybb) *Ubb)/6
      BB_role <- cbind(BB_role,M123_3)

      M123_4 <- (A_2BW*Vbb*(Vbb-1)*Ybb/2 + (C_2BW *(Vbb-1)*Ybb + B_CW * Vbb * (Vbb-1)/2) *Ubb)/7
      BB_role <- cbind(BB_role,rowSums(M123_4),colSums(M123_4))

      M123_5 <- rowSums(A_2BW*Vbb*(Vbb-1)*(Vbb-2)/6 + (C_2BW *(Vbb-1)*(Vbb-2)/2 ) *Ubb)/8 -  M113
      BB_role <- cbind(BB_role,M123_5)


      M221_1 <- rowSums((B_AW *XT + B_ATW *Xbb) *Vbb + C_2BW*Xbb*XT  )/4
      BB_role <- cbind(BB_role,M221_1)

      M221_2 <- ((B_AW *Ubb + A_2BW *Xbb) *Vbb + C_2BW*Xbb*Ubb )/5
      BB_role <- cbind(BB_role,rowSums(M221_2),colSums(M221_2))

      M221_3 <- rowSums(A_2BW*(Ubb-1)*Vbb + C_2BW *Ubb*(Ubb-1)/2)/6 - M211
      BB_role <- cbind(BB_role,M221_3)


      M222_1 <- ((B_AW * Ubb + A_2BW * Xbb) *Ybb*YT +(B_CW*YT+B_CTW*Ybb) *Xbb*Ubb)/5
      BB_role <- cbind(BB_role,rowSums(M222_1),colSums(M222_1))

      M222_2 <- rowSums(A_2BW *(Ubb-1) *Ybb*YT +(B_CW*YT+B_CTW*Ybb) *Ubb*(Ubb-1)/2)/6
      BB_role <- cbind(BB_role,M222_2)

      M222_3 <- ((B_AW * XT + B_ATW * Xbb) *Ybb*Vbb +(B_CW*Vbb + C_2BW*Ybb) *Xbb*XT)/5
      BB_role <- cbind(BB_role,rowSums(M222_3),colSums(M222_3))

      M222_4 <- ((B_AW * Ubb + A_2BW * Xbb) *Ybb*Vbb +(B_CW*Vbb + C_2BW*Ybb) *Xbb*Ubb)/6
      BB_role <- cbind(BB_role,rowSums(M222_4),colSums(M222_4))

      M222_5 <- ((A_2BW * (Ubb-1)*Vbb*YT +(C_2BW*YT + B_CTW* Vbb) *Ubb*(Ubb-1)/2))/7
      BB_role <- cbind(BB_role,rowSums(M222_5),colSums(M222_5))

      M222_6 <- rowSums((B_AW * XT + B_ATW * Xbb) * Vbb*(Vbb-1)/2 + C_2BW* (Vbb-1) *Xbb*XT)/6
      BB_role <- cbind(BB_role,M222_6)

      M222_7 <- ((B_AW * Ubb + A_2BW * Xbb) *Vbb*(Vbb-1)/2 + C_2BW*(Vbb-1) *Xbb*Ubb)/7
      BB_role <- cbind(BB_role,rowSums(M222_7),colSums(M222_7))

      M222_8 <- rowSums(A_2BW *(Ubb-1) *Vbb*(Vbb-1)/2 + C_2BW*(Vbb-1) *Ubb*(Ubb-1)/2)/8 - M212
      BB_role <- cbind(BB_role,M222_8)

      M222_9 <- ((B_AW * Ubb + A_2BW * Xbb) *Vbb*YT +(C_2BW*YT + B_CTW*Vbb) *Xbb*Ubb)/6
      BB_role <- cbind(BB_role,rowSums(M222_9),colSums(M222_9))

      M321_1 <- ((B_AW *(Xbb-1)*XT + B_ATW *Xbb*(Xbb-1)/2) *Vbb + C_2BW *Xbb*(Xbb-1)*XT/2)/5
      BB_role <- cbind(BB_role,rowSums(M321_1),colSums(M321_1))

      M321_2 <- ((A_2BW *XT*(XT-1)/2 + B_ATW *(XT-1)*Ubb) *Vbb + C_2BW *Ubb*XT*(XT-1)/2)/6
      BB_role <- cbind(BB_role,rowSums(M321_2),colSums(M321_2))

      M321_3 <- rowSums((B_AW *Ubb*XT + A_2BW * Xbb*XT + B_ATW *Xbb*Ubb) *Vbb + C_2BW *Xbb*Ubb*XT)/6
      BB_role <- cbind(BB_role,M321_3)

      M321_4 <- ((A_2BW *(Ubb-1)*XT + B_ATW *Ubb*(Ubb-1)/2) *Vbb + C_2BW *Ubb*(Ubb-1)*XT/2)/7
      BB_role <- cbind(BB_role,rowSums(M321_4),colSums(M321_4))

      M321_5 <- rowSums((A_2BW *(Ubb-1)*(Ubb-2)/2) *Vbb + C_2BW *Ubb*(Ubb-1)*(Ubb-2)/6)/8 - M311
      BB_role <- cbind(BB_role,M321_5)

      M131<-M1321<-M1321_1<-M1322<-M1322_1<-M1323<-M1323_1<-M1324<-rep(0,L1)
      M1324_1<-M1325<-M2311<-M2311_1<-M2312<-M2312_1<-M2313<-M2313_1<-rep(0,L1)
      M2314=M2314_1=M2315=M141<-rep(0,L1)

      P0 <- (PW==0)*1
      Q0 <- (QW==0)*1
      if(L1>4){
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


                  role_value <- (PW_F*Q_FCS + QW_F*P_FCS)/6
                  if(sum(role_value)!=0){
                     M131[c(i,j)]<-M131[c(i,j)]+sum(role_value)
                     M131[-c(i,j)]<-M131[-c(i,j)]+role_value
                  }
                  role_value1 <- (PW_F*Q_F1CS*Q_F2CS + (QW_F1*Q_F2CS +QW_F2*Q_F1CS)*P_FCS)/6
                  if(sum(role_value1)!=0){
                     M1321[c(i,j)]<-M1321[c(i,j)]+sum(role_value1)
                     M1321_1[-c(i,j)]<-M1321_1[-c(i,j)]+role_value1
                  }
                  role_value2 <- (PW_F*Q_FCS*Q_F2CS + (QW_F*Q_F2CS +QW_F2*Q_FCS)*P_FCS)/7
                  if(sum(role_value2)!=0){
                     M1322[c(i,j)]<-M1322[c(i,j)]+sum(role_value2)
                     M1322_1[-c(i,j)]<-M1322_1[-c(i,j)]+role_value2
                  }
                  role_value3 <- (PW_F*Q_F1CS*Q_F3CS + (QW_F1*Q_F3CS +QW_F3*Q_F1CS)*P_FCS)/7
                  if(sum(role_value3)!=0){
                     M1323[i]<-M1323[i]+sum(role_value3)
                     M1323[-c(i,j)]<-M1323[-c(i,j)]+role_value3
                     M1323_1[j]<-M1323_1[j]+sum(role_value3)
                  }
                  role_value4 <- (PW_F*Q_FCS*Q_F1CS + (QW_F*Q_F1CS +QW_F1*Q_FCS)*P_FCS)/8
                  if(sum(role_value4)!=0){
                     M1324[c(i,j)]<-M1324[c(i,j)]+sum(role_value4)
                     M1324_1[-c(i,j)]<-M1324_1[-c(i,j)]+role_value4
                  }
                  role_value5 <- (PW_F*Q_FCS*(Q_FCS-1)/2 + QW_F*(Q_FCS-1)*P_FCS)/9
                  if(sum(role_value5)!=0){
                     M1325[c(i,j)]<-M1325[c(i,j)]+sum(role_value5)
                     M1325[-c(i,j)]<-M1325[-c(i,j)]+role_value5
                  }
                  role_value6 <- ((PW_F1*P_F2CS+PW_F2*P_F1CS)*Q_FCS + QW_F*P_F1CS*P_F2CS)/6
                  if(sum(role_value6)!=0){
                     M2311[c(i,j)]<-M2311[c(i,j)]+sum(role_value6)
                     M2311_1[-c(i,j)]<-M2311_1[-c(i,j)]+role_value6
                  }
                  role_value7 <- ((PW_F*P_F2CS+PW_F2*P_FCS)*Q_FCS + QW_F*P_FCS*P_F2CS)/7
                  if(sum(role_value7)!=0){
                     M2312_1[c(i,j)]<-M2312_1[c(i,j)]+sum(role_value7)
                     M2312[-c(i,j)]<-M2312[-c(i,j)]+(role_value7)
                  }
                  role_value8 <- ((PW_F1*P_F3CS+PW_F3*P_F1CS)*Q_FCS + QW_F*P_F1CS*P_F3CS)/7
                  if(sum(role_value8)!=0){
                     M2313[i]<-M2313[i]+sum(role_value8)
                     M2313[-c(i,j)]<-M2313[-c(i,j)]+role_value8
                     M2313_1[j]<-M2313_1[j]+sum(role_value8)
                  }
                  role_value9 <- ((PW_F*P_F3CS+PW_F3*P_FCS)*Q_FCS + QW_F*P_FCS*P_F3CS)/8
                  if(sum(role_value9)!=0){
                     M2314[i]<-M2314[i]+sum(role_value9)
                     M2314_1[j]<-M2314_1[j]+sum(role_value9)
                     M2314_1[-c(i,j)]<-M2314_1[-c(i,j)]+role_value9
                  }
                  role_value10 <- (PW_F*(P_FCS-1)*Q_FCS + QW_F*P_FCS*(P_FCS-1)/2)/9
                  if(sum(role_value10)!=0){
                     M2315[c(i,j)]<-M2315[c(i,j)]+sum(role_value10)
                     M2315[-c(i,j)]<-M2315[-c(i,j)]+role_value10
                  }
               }
               if(i<j){
                  for(k in (j):L1){
                     if(j<k){
                        P_F4 <- P[,i]*P[,j]*P[,k]*P[,-c(i,j,k)]
                        Q_F4 <- Q[i,]*Q[j,]*Q[k,]*t(Q[-c(i,j,k),])

                        P_F4CS <- colSums(P_F4)
                        Q_F4CS <- colSums(Q_F4)


                        PW_ijkz <- PW[,i]+PW[,j]+PW[,k]+PW[,-c(i,j,k)]
                        QW_ijkz <- QW[i,]+QW[j,]+QW[k,]+t(QW[-c(i,j,k),])

                        PW_F4 <- colSums(PW_ijkz*P_F4)
                        QW_F4 <- colSums(QW_ijkz*Q_F4)

                        role_value11 <- (PW_F4 * Q_F4CS + QW_F4 * P_F4CS)/8##########################
                        if(sum(role_value11)!=0){
                           M141[c(i,j,k)]<-M141[c(i,j,k)]+sum(role_value11)
                           M141[-c(i,j,k)]<-M141[-c(i,j,k)]+role_value11
                        }
                     }
                  }
               }
            }
         }
         M141<-M141/4
      }
      else{
         for(i in 1:4){
            for(j in 1:4){
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


                  role_value <- (PW_F*Q_FCS + QW_F*P_FCS)/6
                  if(sum(role_value)!=0){
                     M131[c(i,j)]<-M131[c(i,j)]+sum(role_value)
                     M131[-c(i,j)]<-M131[-c(i,j)]+role_value
                  }
                  role_value1 <- (PW_F*Q_F1CS*Q_F2CS + (QW_F1*Q_F2CS +QW_F2*Q_F1CS)*P_FCS)/6
                  if(sum(role_value1)!=0){
                     M1321[c(i,j)]<-M1321[c(i,j)]+sum(role_value1)
                     M1321_1[-c(i,j)]<-M1321_1[-c(i,j)]+role_value1
                  }
                  role_value2 <- (PW_F*Q_FCS*Q_F2CS + (QW_F*Q_F2CS +QW_F2*Q_FCS)*P_FCS)/7
                  if(sum(role_value2)!=0){
                     M1322[c(i,j)]<-M1322[c(i,j)]+sum(role_value2)
                     M1322_1[-c(i,j)]<-M1322_1[-c(i,j)]+role_value2
                  }
                  role_value3 <- (PW_F*Q_F1CS*Q_F3CS + (QW_F1*Q_F3CS +QW_F3*Q_F1CS)*P_FCS)/7
                  if(sum(role_value3)!=0){
                     M1323[i]<-M1323[i]+sum(role_value3)
                     M1323[-c(i,j)]<-M1323[-c(i,j)]+role_value3
                     M1323_1[j]<-M1323_1[j]+sum(role_value3)
                  }
                  role_value4 <- (PW_F*Q_FCS*Q_F1CS + (QW_F*Q_F1CS +QW_F1*Q_FCS)*P_FCS)/8
                  if(sum(role_value4)!=0){
                     M1324[c(i,j)]<-M1324[c(i,j)]+sum(role_value4)
                     M1324_1[-c(i,j)]<-M1324_1[-c(i,j)]+role_value4
                  }
                  role_value5 <- (PW_F*Q_FCS*(Q_FCS-1)/2 + QW_F*(Q_FCS-1)*P_FCS)/9
                  if(sum(role_value5)!=0){
                     M1325[c(i,j)]<-M1325[c(i,j)]+sum(role_value5)
                     M1325[-c(i,j)]<-M1325[-c(i,j)]+role_value5
                  }
                  role_value6 <- ((PW_F1*P_F2CS+PW_F2*P_F1CS)*Q_FCS + QW_F*P_F1CS*P_F2CS)/6
                  if(sum(role_value6)!=0){
                     M2311[c(i,j)]<-M2311[c(i,j)]+sum(role_value6)
                     M2311_1[-c(i,j)]<-M2311_1[-c(i,j)]+role_value6
                  }
                  role_value7 <- ((PW_F*P_F2CS+PW_F2*P_FCS)*Q_FCS + QW_F*P_FCS*P_F2CS)/7
                  if(sum(role_value7)!=0){
                     M2312_1[c(i,j)]<-M2312_1[c(i,j)]+sum(role_value7)
                     M2312[-c(i,j)]<-M2312[-c(i,j)]+(role_value7)
                  }
                  role_value8 <- ((PW_F1*P_F3CS+PW_F3*P_F1CS)*Q_FCS + QW_F*P_F1CS*P_F3CS)/7
                  if(sum(role_value8)!=0){
                     M2313[i]<-M2313[i]+sum(role_value8)
                     M2313[-c(i,j)]<-M2313[-c(i,j)]+role_value8
                     M2313_1[j]<-M2313_1[j]+sum(role_value8)
                  }
                  role_value9 <- ((PW_F*P_F3CS+PW_F3*P_FCS)*Q_FCS + QW_F*P_FCS*P_F3CS)/8
                  if(sum(role_value9)!=0){
                     M2314[i]<-M2314[i]+sum(role_value9)
                     M2314_1[j]<-M2314_1[j]+sum(role_value9)
                     M2314_1[-c(i,j)]<-M2314_1[-c(i,j)]+role_value9
                  }
                  role_value10 <- (PW_F*(P_FCS-1)*Q_FCS + QW_F*P_FCS*(P_FCS-1)/2)/9
                  if(sum(role_value10)!=0){
                     M2315[c(i,j)]<-M2315[c(i,j)]+sum(role_value10)
                     M2315[-c(i,j)]<-M2315[-c(i,j)]+role_value10
                  }
               }
               M141 <- (sum(PW[,1]+PW[,2]+PW[,3]+PW[,4])*sum(Q[1,]*Q[2,]*Q[3,]*Q[4,]) + sum(QW[1,]+QW[2,]+QW[3,]+Q[4,])*sum(P[,1]*P[,2]*PW[,3]*P[,4]))/8
            }
         }
      }
      BB_role<-cbind(BB_role,M131/6,M1321/2,M1321_1/2,M1322_1/2,M1322/2,M1323/2,M1323_1/2,M1324/2,M1324_1/2,M1325/6,M2311/2,M2311_1/2,M2312/2,M2312_1/2,M2313/2,M2313_1/2,M2314/2,M2314_1/2,M2315/6,M141)
      colnames(BB_role)<-paste0("role",c(1:70))
      weighted_role <- role
      weighted_role[rownames(BB_role),]<-BB_role
      colnames(weighted_role)<-paste0("weighted role",c(1:70))

      mean_role <- weighted_role/role
      mean_role <- replace(mean_role, which(is.nan(mean_role)), NA)

      return(mean_role)
   }
   else
      return(role)
}

