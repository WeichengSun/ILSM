#' Counts of interconnection motifs
#'
#' Counting the frequencies of interconnection motifs for a tripartite interaction network.
#'
#' @param network.or.subnet_mat1 Igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork.
#' @param subnet_mat2 The matrix representing one subnetwork.
#' @param weighted Logical. Default to FALSE. If TRUE, the counts of motifs are weighted by link strengths in each motif.

#' @import igraph
#'
#' @export
#'
#' @details
#' In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.
#' An interconnection motif is defined to comprise three sets of connected nodes: the connector nodes (belonging to b-nodes), the nodes in one subnetwork (belonging to a-nodes in the P subnetwork), and the nodes in the other subnetwork (belonging to c-nodes in the Q subnetwork). Each motif has maximumly 6 nodes, resulting in a total of 48 distinct motif forms.
#' The algorithm for counting interconnection motifs is designed by extending the fast approach from Simmons et al.(2019), which uses mathematical operations directly on the bi-adjacency matrix.
#' \strong{weighted}
#'
#' If the \code{weighted} = FALSE, the input network can be an "igraph" object or two matrices. If a weighted network is provided, it will be transformed to a binary network.
#' If the \code{weighted} = TRUE, the input network can only be two matrices. Weighted motifs are considered, then the counts are weighed by the link strength in each motif.
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
#'  return a numeric vector of the number of 48 interconnection motifs.  The motifs are named ???MABC-i???: M means ???motif???, ???A??? is the number of a-nodes, ???B??? is the number of b-nodes, ???C??? is the number of c-nodes and ???i??? is the serial number for the motifs with the same ???ABC???. See 'Multi_motif' for the forms.
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
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#'#'
#' @examples
#'
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#' icmotif_count(Net)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' icmotif_count(Net)
#'
#'
#' set.seed(12)
#' MAT <- build_net(11,22,21,0.2,output_matrices=TRUE)
#' icmotif_count(MAT[["subnet_P"]],MAT[["subnet_Q"]])
#'
#'##input as binary matrices,with row names.
#' md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
#' md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
#' icmotif_count(md1,md2)
#'
#'
#'##input as weighted matrices,with row numbers as row names.
#' mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
#' mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
#' icmotif_count(mdw1,mdw2,weighted=T)
#'
icmotif_count <-function(network.or.subnet_mat1, subnet_mat2=NULL, weighted = FALSE){
   if(weighted){
      if(inherits(network.or.subnet_mat1,"igraph")==T){
         network<-adjust_net(network.or.subnet_mat1)
         PHP<-as.matrix(network[])
         dimnames(PHP)<-NULL
         PH<-PHP[(V(network)$level)==0,(V(network)$level)==1]
         HP<-PHP[(V(network)$level)==1,(V(network)$level)==2]
      }
      else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
         mat1<-network.or.subnet_mat1
         mat2<-subnet_mat2
         if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
            rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
            rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         }
         if(nrow(mat1)!=nrow(mat2))
            message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
         if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
            matrow<-unique(c(rownames(mat1),rownames(mat2)))
         else
            stop("Make sure matrices either have no row names or have full row names. No NA!!!")
         mat_1<-matrix(0,length(matrow),ncol(mat1))
         rownames(mat_1)<-matrow
         mat_1[rownames(mat1),]<-mat1
         mat_2<-matrix(0,length(matrow),ncol(mat2))
         rownames(mat_2)<-matrow
         mat_2[rownames(mat2),]<-mat2
         mat1<-mat_1
         mat2<-mat_2
         dimnames(mat1)<-NULL
         dimnames(mat2)<-NULL
         PH<-t(mat1)
         HP<-mat2
         logi<- colSums(PH) * rowSums(HP) != 0
         #(apply(PH,2,sum)*apply(HP,1,sum))!=0
         PH<-PH[,logi]
         HP<-HP[logi,]
      }
      else
         stop("Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
      L1<-ncol(PH)
      if(L1 < 4L)
         stop("Error: please input a large 'number of interconnecting species >=4' network data!!!")
      Ob<-colSums(PH)
      Rb<-rowSums(HP)
      Ubb<-t(PH)%*%PH
      Xbb<-t(PH)%*%(PH==0)
      XT<-t(Xbb)
      Vbb<-HP%*%t(HP)
      Ybb<-HP%*%(t(HP)==0)
      YT<-t(Ybb)
      UVbb<-Ubb*Vbb
      Two<-function(x){ ((sum(x))^2 - sum(x^2)) / 2 }
      Three<-function(x){ s1 <- sum(x); s2 <- sum(x^2); s3 <- sum(x^3)
      ((s1^3) - 3*s1*s2 + 2*s3) / 6}
      Four<-function(x){ s1 <- sum(x); s2 <- sum(x^2); s3 <- sum(x^3);
      s4 <- sum(x^4)
      ((s1^4) - 6*(s1^2)*s2 + 3*(s2^2) + 8*s1*s3 - 6*s4)/ 24}
      mat_tri<-function(mat){ sum( mat[upper.tri(mat)] ) }

      PP<-PH%*%HP

      M111 <- sum(PP) # sum(Ob*Rb)

      M112 <- sum(Ob*apply(HP,1,Two)) #=sum(PH%*%apply(HP,1,Two))

      M113 <- sum(Ob*apply(HP,1,Three))  #=sum(PH%*%apply(HP,1,Two_3))

      M114 <- sum(Ob*apply(HP,1,Four)) #=sum(PH%*%apply(HP,1,function(s){sum(s)*(sum(s)-1)*(sum(s)-2)*(sum(s)-3)/24}))

      M211 <- sum(apply((PH),2,Two)*Rb) #=sum(apply((PH),2,Two)%*%HP)

      M212 <- sum(apply((PH),2,Two)*apply(HP,1,Two))  #=sum(apply(PH,2,Two)%*%apply(HP,1,Two))

      M213 <- sum(apply((PH),2,Two)*apply(HP,1,Three))  #=sum(apply(PH,2,Two)%*%apply(HP, 1,Two_3))

      M311 <- sum(apply((PH),2,Three)*Rb)  #=sum(t(HP)%*%apply(t(PH),1,Two_3))

      M312 <- sum(apply((PH),2,Three)*apply(HP,1,Two))  #=sum(apply(PH,2,Two_3)%*%apply(HP,1, Two))

      M411 <- sum(apply((PH),2,Four)*Rb) #=sum(t(HP)%*%apply(t(PH),1,function(s){sum(s)*(sum(s)-1)*(sum(s)-2)*(sum(s)-3)/24}))

      M121 <- mat_tri(UVbb)   #(sum(UVbb)-M111)/2

      M122_1 <- sum(Ubb*Ybb*YT)/2 #=sum(PH_add*HP_ratio*t(HP_ratio))/2

      M122_2 <- sum(UVbb*Ybb) #=sum(PHP_add*HP_ratio)##

      Two_Ubb <- (Ubb^2 - t(apply(PH,2,function(x){colSums((x*PH)^2)})))/2
      Two_Vbb <- (Vbb^2 - t(apply(HP,1,function(x){colSums((x*t(HP))^2)})))/2

      M122_3 <- mat_tri(Ubb*Two_Vbb)  #(sum(Ubb*Two(Vbb))-M112)/2

      M123_1 <- sum(Ubb*YT*(Ybb^2 - t(apply(HP,1,function(x){colSums((x*t(HP==0))^2)})))/2)  #sum(Ubb*Two(Ybb)*YT)

      M123_2 <- sum(UVbb*(YT^2-apply(HP,1,function(x){colSums((x*t(HP==0))^2)}))/2)  #sum(UVbb*Two(YT))

      M123_3 <- sum(UVbb*Ybb*YT)/2

      M123_4 <- sum(Ubb*Two_Vbb*Ybb) #sum(Ubb*Two(Vbb)*Ybb)

      M123_5 <- mat_tri(Ubb*(Vbb^3-3*Vbb*t(apply(HP,1,function(x){colSums((x*t(HP))^2)}))+2*t(apply(HP,1,function(x){colSums((x*t(HP))^3)})))/6)  #(sum(Ubb*Three(Vbb))-M113)/2

      M221_1 <- sum(Xbb*XT*Vbb)/2 #mat_tri(Xbb*XT*Vbb)

      M221_2 <- sum(UVbb*Xbb)

      M221_3 <- mat_tri(Two_Ubb *Vbb)  ###(sum(Two(Ubb)*Vbb)-M211)/2

      M222_1 <- sum(Ubb*Xbb*Ybb*YT)

      M222_2 <- sum(Two_Ubb*Ybb*YT)/2  ##sum(Two(Ubb)*Ybb*YT)/2

      ###########################################


      M222_3 <- sum(Xbb*XT*Ybb*Vbb)



      M222_4 <- sum(UVbb*Xbb*Ybb)


      M222_5 <- sum(Two_Ubb*Vbb*YT)    #sum(Two(Ubb)*Vbb*YT)


      M222_6 <- sum(Xbb*XT*Two_Vbb)/2 # sum(Xbb*XT*Two(Vbb))/2


      M222_7 <- sum(Ubb*Xbb*Two_Vbb) #sum(Ubb*Xbb*Two(Vbb))



      M222_8 <- mat_tri(Two_Ubb * Two_Ubb) # (sum(Two(Ubb)*Two(Vbb))-M212)/2



      M222_9 <- sum(UVbb*Xbb*YT)



      M321_1 <- sum((Xbb^2 - t(apply(PH,2,function(x){colSums((x*(PH==0))^2)})))/2 *XT*Vbb)   # sum(Two(Xbb)*XT*Vbb)


      M321_2 <- sum(UVbb*(XT^2-apply(PH,2,function(x){colSums((x*(PH==0))^2)}))/2)  # sum(UVbb*Two(XT))



      M321_3 <- sum(UVbb*Xbb*XT)/2


      M321_4 <- sum(Two_Ubb*XT*Vbb)   # sum(Two(Ubb)*XT*Vbb)



      M321_5 <- mat_tri((Ubb^3-3*Ubb*t(apply(PH,2,function(x){colSums((x*PH)^2)}))+2*t(apply(PH,2,function(x){colSums((x*PH)^3)})))/6*Vbb) #(sum(Three(Ubb)*Vbb)-M311)/2


      M131 <-
         sum((PP^3 -
                 3 * PP * t(apply(PH, 1, function(x){colSums((x * HP) ^ 2) })) +
                 2 * t(apply(PH, 1, function(x) {colSums((x * HP) ^ 3)})))
             / 6)

      # apply(PH, 1, function(x){{
      #    colSums((x*HP)^2)
      # }})

      M141 <-
         sum((PP^4 +
                 8 * PP * t(apply(PH, 1, function(x) {colSums((x * HP) ^ 3)})) -
                 6 * PP^2 * t(apply(PH, 1, function(x) {colSums((x * HP) ^ 2)}))-
                 3 * t(apply(PH, 1, function(x) {colSums((x * HP) ^ 4)}))+
                 6 * t(apply(PH, 1, function(x) {apply((x * HP) ^ 2,2,Two)}))
         ) / 24)


      M1321<-M1322<-M1323<-M1324<-M1325<-M2311<-M2312<-M2313<-M2314<-M2315<-0
      for(i in 1:L1){
         for(j in 1:L1){
            if(i!=j){
               PH_F<-colSums(PH[,i]*PH[,j]*(PH[,-c(i,j)]))
               PH_F_2T <- colSums((PH[,i]*PH[,j]*(PH[,-c(i,j)]))^2)
               PH_F1<-colSums(PH[,i]*PH[,j]*(PH[,-c(i,j)]==0))
               PH_F2<-colSums((PH[,i]==0)*(PH[,j]==0)*(PH[,-c(i,j)]))
               PH_F3<-colSums((PH[,i]==0)*(PH[,j])*(PH[,-c(i,j)]))
               HP_F<-colSums(HP[i,]*HP[j,]*t(HP[-c(i,j),]))
               HP_F_2T <- colSums((HP[i,]*HP[j,]*t(HP[-c(i,j),]))^2)
               HP_F1<-colSums(HP[i,]*HP[j,]*t(HP[-c(i,j),]==0))
               HP_F2<-colSums((HP[i,]==0)*(HP[j,]==0)*t(HP[-c(i,j),]))
               HP_F3<-colSums((HP[i,]==0)*(HP[j,])*t(HP[-c(i,j),]))
               ########
               M1321<-M1321+sum(PH_F*HP_F1*HP_F2)
               M1322<-M1322+sum(PH_F*HP_F*HP_F2)
               M1323<-M1323+sum(PH_F*HP_F1*HP_F3)
               M1324<-M1324+sum(PH_F*HP_F*HP_F1)
               M1325<-M1325+sum(PH_F*((HP_F^2-HP_F_2T)/2))
               M2311<-M2311+sum(PH_F1*PH_F2*HP_F)
               M2312<-M2312+sum(PH_F*PH_F2*HP_F)
               M2313<-M2313+sum(PH_F1*PH_F3*HP_F)
               M2314<-M2314+sum(PH_F*PH_F3*HP_F)
               M2315<-M2315+sum(((PH_F^2-PH_F_2T)/2)*HP_F)
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
   }
   else
      ###################################################################################
   {

   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network<-adjust_net(network.or.subnet_mat1)
      PHP<-as.matrix(network[])
      PHP[PHP>0]<-1
      dimnames(PHP)<-NULL
      PH<-PHP[(V(network)$level)==0,(V(network)$level)==1]
      HP<-PHP[(V(network)$level)==1,(V(network)$level)==2]
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
      mat1<-network.or.subnet_mat1
      mat1[mat1>0]<-1
      mat2<-subnet_mat2
      mat2[mat2>0]<-1
      if(is.null(rownames(mat1)) | is.null(rownames(mat2))){
         rownames(mat1)<-paste0("mid_spe",seq=1:nrow(mat1))
         rownames(mat2)<-paste0("mid_spe",seq=1:nrow(mat2))
         matrow<-unique(c(rownames(mat1),rownames(mat2)))
      }
      if(nrow(mat1)!=nrow(mat2))
         message("re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
      if(!is.null(rownames(mat1)) & !is.null(rownames(mat2)) & sum(is.na(rownames(mat1)))==0 & sum(is.na(rownames(mat2)))==0)
         matrow<-unique(c(rownames(mat1),rownames(mat2)))
      else
         stop("Make sure matrices either have no row names or have full row names. No NA!!!")
      mat_1<-matrix(0,length(matrow),ncol(mat1))
      rownames(mat_1)<-matrow
      mat_1[rownames(mat1),]<-mat1
      mat_1[mat_1>0]<-1
      mat_2<-matrix(0,length(matrow),ncol(mat2))
      rownames(mat_2)<-matrow
      mat_2[rownames(mat2),]<-mat2
      mat_2[mat_2>0]<-1
      mat1<-mat_1
      mat2<-mat_2
      dimnames(mat1)<-NULL
      dimnames(mat2)<-NULL
      PH<-t(mat1)
      HP<-mat2
      logi<-(apply(PH,2,sum)*apply(HP,1,sum))!=0
      PH<-PH[,logi]
      HP<-HP[logi,]
   }
   else
      stop("Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   L1<-ncol(PH)
   if(L1 < 4L)
      stop("Error: please input a large 'number of interconnecting species >=4' network data!!!")
   Ob<-colSums(PH)
   Rb<-rowSums(HP)
   Ubb<-t(PH)%*%PH
   Xbb<-t(PH)%*%(1-PH)
   XT<-t(Xbb)
   Vbb<-HP%*%t(HP)
   Ybb<-HP%*%(1-t(HP))
   YT<-t(Ybb)
   UVbb<-Ubb*Vbb
   Two<-function(a){ return(a*(a-1)/2) }
   Three<-function(a){ return(a*(a-1)*(a-2)/6)}
   Four<-function(a){ return(a*(a-1)*(a-2)*(a-3)/24)}

   PP<-PH%*%HP

   M111 <- sum(Ob*Rb) #=sum(PP)

   M112 <- sum(Ob*Two(Rb)) #=sum(PH%*%apply(HP,1,Two))

   M113 <- sum(Ob*Three(Rb))  #=sum(PH%*%apply(HP,1,Two_3))

   M114 <- sum(Ob*Four(Rb)) #=sum(PH%*%apply(HP,1,function(s){sum(s)*(sum(s)-1)*(sum(s)-2)*(sum(s)-3)/24}))

   M211 <- sum(Two(Ob)*Rb) #=sum(apply((PH),2,Two)%*%HP)

   M212 <- sum(Two(Ob)*Two(Rb))  #=sum(apply(PH,2,Two)%*%apply(HP,1,Two))

   M213 <- sum(Two(Ob)*Three(Rb))  #=sum(apply(PH,2,Two)%*%apply(HP, 1,Two_3))

   M311 <- sum(Three(Ob)*Rb)  #=sum(t(HP)%*%apply(t(PH),1,Two_3))

   M312 <- sum(Three(Ob)*Two(Rb))  #=sum(apply(PH,2,Two_3)%*%apply(HP,1, Two))

   M411 <- sum(Four(Ob)*Rb) #=sum(t(HP)%*%apply(t(PH),1,function(s){sum(s)*(sum(s)-1)*(sum(s)-2)*(sum(s)-3)/24}))

   M121 <- (sum(UVbb)-M111)/2 #=sum(PP*(PP-1)/2)##

   M122_1 <- sum(Ubb*Ybb*YT)/2 #=sum(PH_add*HP_ratio*t(HP_ratio))/2

   M122_2 <- sum(UVbb*Ybb) #=sum(PHP_add*HP_ratio)##

   M122_3 <- (sum(Ubb*Two(Vbb))-M112)/2 # M<-PHP_add*(HP_add-1); =sum(M[lower.tri(M)])/2

   M123_1 <- sum(Ubb*Two(Ybb)*YT) #=sum(PH_add*(HP_ratio*(HP_ratio-1)/2)*t(HP_ratio))

   M123_2 <- sum(UVbb*Two(YT)) #=sum(PHP_add*(HP_ratio*(HP_ratio-1)/2))##

   M123_3 <- sum(UVbb*Ybb*YT)/2 #=sum(PHP_add*HP_ratio*t(HP_ratio))/2##

   M123_4 <- sum(Ubb*Two(Vbb)*Ybb) #=sum(PHP_add*((HP_add-1)/2)*HP_ratio)##

   M123_5 <- (sum(Ubb*Three(Vbb))-M113)/2 #  M<-PHP_add*((HP_add-1)*(HP_add-2)/6); =sum(M[lower.tri(M)])

   M221_1 <- sum(Xbb*XT*Vbb)/2 #=sum(HP_add*PH_ratio*t(PH_ratio))/2

   M221_2 <- sum(UVbb*Xbb) #=sum(PHP_add*(PH_ratio))##

   M221_3 <- (sum(Two(Ubb)*Vbb)-M211)/2 #M<-PHP_add*(PH_add-1); =sum(M[lower.tri(M)])/2

   M222_1 <- sum(Ubb*Xbb*Ybb*YT) #=sum(PH_add*PH_ratio*HP_ratio*t(HP_ratio))

   M222_2 <- sum(Two(Ubb)*Ybb*YT)/2 #=sum(PH_add*(PH_add-1)*HP_ratio*t(HP_ratio))/4

   M222_3 <- sum(Xbb*XT*Ybb*Vbb) #=sum(PH_ratio*t(PH_ratio)*HP_add*HP_ratio)

   M222_4 <- sum(UVbb*Xbb*Ybb) #=sum(PHP_add*PH_ratio*HP_ratio)##

   M222_5 <- sum(Two(Ubb)*Vbb*YT) #=sum((PHP_add*(PH_add-1)/2)*HP_ratio)##

   M222_6 <- sum(Xbb*XT*Two(Vbb))/2 #=sum(PH_ratio*t(PH_ratio)*HP_add*(HP_add-1))/4

   M222_7 <- sum(Ubb*Xbb*Two(Vbb)) #=sum(PHP_add*PH_ratio*(HP_add-1)/2)##

   M222_8 <- (sum(Two(Ubb)*Two(Vbb))-M212)/2 #M<-PHP_add*(PH_add-1)*(HP_add-1)/4; =sum(M[lower.tri(M)])

   M222_9 <- sum(UVbb*Xbb*YT) #=sum(PHP_add*PH_ratio*t(HP_ratio))##

   M321_1 <- sum(Two(Xbb)*XT*Vbb) #=sum((PH_ratio*(PH_ratio-1)/2)*t(PH_ratio)*HP_add)

   M321_2 <- sum(UVbb*Two(XT)) #=sum((PH_ratio*(PH_ratio-1)/2)*PHP_add)##

   M321_3 <- sum(UVbb*Xbb*XT)/2 #=sum(PHP_add*PH_ratio*t(PH_ratio))/2##

   M321_4 <- sum(Two(Ubb)*XT*Vbb) #=sum((PHP_add*(PH_add-1)/2)*PH_ratio)##

   M321_5 <- (sum(Three(Ubb)*Vbb)-M311)/2 #M<-(PHP_add*(PH_add-1)*(PH_add-2)/6); =sum(M[lower.tri(M)])


   M131 <- sum(Three(PP))

   M141 <- sum(Four(PP))

   M1321<-M1322<-M1323<-M1324<-M1325<-M2311<-M2312<-M2313<-M2314<-M2315<-0
   for(i in 1:L1){
      for(j in 1:L1){
         if(i!=j){
            PH_F<-colSums(PH[,i]*PH[,j]*(PH[,-c(i,j)]))
            PH_F1<-colSums(PH[,i]*PH[,j]*(1-PH[,-c(i,j)]))
            PH_F2<-colSums((1-PH[,i])*(1-PH[,j])*(PH[,-c(i,j)]))
            PH_F3<-colSums((1-PH[,i])*(PH[,j])*(PH[,-c(i,j)]))
            HP_F<-colSums(HP[i,]*HP[j,]*t(HP[-c(i,j),]))
            HP_F1<-colSums(HP[i,]*HP[j,]*t(1-HP[-c(i,j),]))
            HP_F2<-colSums((1-HP[i,])*(1-HP[j,])*t(HP[-c(i,j),]))
            HP_F3<-colSums((1-HP[i,])*(HP[j,])*t(HP[-c(i,j),]))
            ########
            M1321<-M1321+sum(PH_F*HP_F1*HP_F2)
            M1322<-M1322+sum(PH_F*HP_F*HP_F2)
            M1323<-M1323+sum(PH_F*HP_F1*HP_F3)
            M1324<-M1324+sum(PH_F*HP_F*HP_F1)
            M1325<-M1325+sum(PH_F*HP_F*(HP_F-1)/2)
            M2311<-M2311+sum(PH_F1*PH_F2*HP_F)
            M2312<-M2312+sum(PH_F*PH_F2*HP_F)
            M2313<-M2313+sum(PH_F1*PH_F3*HP_F)
            M2314<-M2314+sum(PH_F*PH_F3*HP_F)
            M2315<-M2315+sum((PH_F*(PH_F-1)/2)*HP_F)
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
   }
   motif<-c(M111, M112, M113, M114, M211, M212, M213, M311, M312, M411,
            M121, M122_1, M122_2, M122_3, M123_1, M123_2, M123_3,
            M123_4, M123_5, M221_1, M221_2, M221_3, M222_1, M222_2,
            M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_9, M321_1,
            M321_2, M321_3, M321_4, M321_5, M131, M1321, M1322,
            M1323, M1324, M1325, M2311, M2312, M2313, M2314, M2315, M141)
   names(motif)<-c("M111", "M112","M113", "M114", "M211", "M212", "M213", "M311", "M312", "M411",
     "M121", "M122-1", "M122-2", "M122-3", "M123-1", "M123-2", "M123-3",
     "M123-4", "M123-5", "M221-1", "M221-2", "M221-3", "M222-1", "M222-2",
     "M222-3", "M222-4", "M222-5", "M222-6", "M222-7", "M222-8", "M222-9", "M321-1",
     "M321-2", "M321-3", "M321-4", "M321-5", "M131", "M132-1", "M132-2",
     "M132-3", "M132-4", "M132-5", "M231-1", "M231-2", "M231-3", "M231-4", "M231-5", "M141")
   return(motif)
}
