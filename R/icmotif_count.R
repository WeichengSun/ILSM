#' Calculating the number of 48 motifs
#'
#' Calculating the number of 48 motifs from a tripartite interaction network.
#'
#' @param network.or.subnet_mat1 Either a tripartite network of 'igraph' class which contains three groups of species and interactions within subnetwork without interactions between each group of species, or a numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#' @param ic_number The max number of interconnection code in each motif. Defaults to 4 because of containing all motifs, also can input 2 to calculate the top 36 for rapid computing speed.
#'
#' @import igraph
#'
#' @export
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
#' If the type of inputting is data frame or matrix, please make sure the row of \code{network.or.subnet_mat1} and \code{subnet_mat2} correspond with the second group of species that both belong to two subnetworks and interact with other groups of species.
#' \itemize{
#' \item Try to make the rows of both matrices have the same attributes. Or we default:
#'
#' \item When the two matrices can have different numbers of rows:
#' \itemize{
#' \item 1. If both matrices have row names, then the function counts all row names to produce two new matrices with the same row names.
#' \item 2. If at most one matrix has row names, the function assigns new row names to both matrices on a row-to-row basis (any extra row names are assigned a new value) and then counts all row names to produce two new matrices with the same row names.
#' }
#'
#' \item When the two matrices can have the same numbers of rows:
#' \itemize{
#' \item No matter how the row names of the two matrices are arranged, as long as the row names are exactly the same; But we don't handle matrices with empty row names (the function will give an error).
#' }
#'
#' \item The two matrices can have different numbers of rows, but read our default handling carefully to make sure the calculation is accurate when using this function!!!
#' }
#' About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}.
#'
#'
#'
#' @return
#' If \code{ic_number} = 4,
#'  return a numeric vector with the number of 48 motifs: M111, M112, M113,
#'  M114, M211, M212, M213, M311, M312, M411, M121_1, M122_1, M122_2, M122_3,
#'  M123_1, M123_2, M123_3, M123_4, M123_5, M221_1, M221_2, M221_3, M222_1,
#'  M222_2, M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_8, M321_1,
#'  M321_2, M321_3, M321_4, M321_5, M131, M132-1, M132-2, M132-3, M132-4,
#'  M132-5, M231-1, M231-2, M231-3, M231-4, M231-5, M141.
#'
#' However if \code{ic_nmeber} = 2,
#'  only return a numeric vector with the number of top 36 motifs: M111, M112,
#'  M113, M114, M211, M212, M213, M311, M312, M411, M121_1, M122_1, M122_2,
#'  M122_3, M123_1, M123_2, M123_3, M123_4, M123_5, M221_1, M221_2, M221_3,
#'  M222_1, M222_2, M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_8,
#'  M321_1, M321_2, M321_3, M321_4, M321_5
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
#' Pilosof, S., Porter, M. A., Pascual, M., & Kéfi, S. (2017). The multilayer nature of ecological networks. Nature Ecology & Evolution, 1(4), 0101.
#'
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#'
#'
#' @examples
#'
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' icmotif_count(Net)
#' icmotif_count(Net, ic_number=2)
#'
#' set.seed(12)
#' MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
#' icmotif_count(MAT[[3]],MAT[[4]])
#'
#' md1 <- matrix(sample(c(0,1),120,replace=TRUE),8,15)
#' md2 <- matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' icmotif_count(md1,md2)
#'
#' R <- rownames(MAT[[4]])[12]
#' MR <- MAT[[4]][12,]
#' MAT[[4]] <- MAT[[4]][-12,]
#' MAT[[4]] <- rbind(MAT[[4]],MR)
#' rownames(MAT[[4]])[22] <- R
#'
#' icmotif_count(MAT[[3]],MAT[[4]])
#'
#'

icmotif_count <- function(network.or.subnet_mat1, subnet_mat2=NULL, ic_number=4){
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

   if(ic_number==2)
      return(c(M111, M112, M113, M114, M211, M212, M213, M311, M312, M411,
               M121, M122_1, M122_2, M122_3, M123_1, M123_2, M123_3,
               M123_4, M123_5, M221_1, M221_2, M221_3, M222_1, M222_2,
               M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_9, M321_1,
               M321_2, M321_3, M321_4, M321_5))

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
   motif<-c(M111, M112, M113, M114, M211, M212, M213, M311, M312, M411,
            M121, M122_1, M122_2, M122_3, M123_1, M123_2, M123_3,
            M123_4, M123_5, M221_1, M221_2, M221_3, M222_1, M222_2,
            M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M222_9, M321_1,
            M321_2, M321_3, M321_4, M321_5, M131, M1321, M1322,
            M1323, M1324, M1325, M2311, M2312, M2313, M2314, M2315, M141)
   ###########################################################
   # if(subnet_motif){
   #    subnet1_motif<-bmotif::mcount(PH,six_node = T, normalisation = TRUE, mean_weight = F, standard_dev = F)
   #    subnet2_motif<-bmotif::mcount(HP,six_node = T, normalisation = TRUE, mean_weight = F, standard_dev = F)
   #    subnet1_node_position<-bmotif::node_positions(PH,six_node=TRUE,weights_method="none")
   #    subnet2_node_position<-bmotif::node_positions(HP,six_node=TRUE,weights_method="none")
   #    subnet1_link_position<-bmotif::link_positions(PH,six_node=TRUE,weights =FALSE)
   #    subnet2_link_position<-bmotif::link_positions(HP,six_node=TRUE,weights =FALSE)
   #    motif_list<-list()
   #    motif_list$tripartitenet_motif<-motif
   #    motif_list$subnetwork1_motif<-subnet1_motif
   #    motif_list$sunnetwork2_motif<-subnet2_motif
   #    motif_list$subnet1_node_position<-subnet1_node_position
   #    motif_list$subnet2_node_position<-subnet2_node_position
   #    motif_list$subnet1_link_position<-subnet1_link_position
   #    motif_list$subnet2_link_position<-subnet2_link_position
   #
   #    return(motif_list)
   # }
   # else
      return(motif)
}
