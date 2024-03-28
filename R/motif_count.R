#' Calculating the number of 44 motifs
#'
#' Calculating the number of 44 motifs from a tripartite interaction network.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class, or a numeric matrix(or data.frame) representing interactions between two groups of species. The network contains
#'  interlayer links and without intralayer links. Each row and column of matrix represents the species in the first and second layers of the tripartite network respectively.
#'  Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents the species in the second and third layers of the
#'  tripartite network respectively. Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#'
#' @param subnet_motif Logical. Whether to calculate the structure indices of the two bipartite networks dismantled by this multilayer network based on the principle of the package bmotif. Default to FALSE: NO.
#' @import igraph
#'
#' @export
#' @details
#'
#' \strong{network.or.subnet_mat1}
#'
#' There are two types of data that can be processed:
#' \itemize{
#' \item{Input in a network of type "igraph" alone.}
#'
#' \item{Must be entered as data frame or matrix with \code{subnet_mat2}}}
#'
#' About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}
#'
#'
#'
#' @return
#' If \code{subnet_motif} = FALSE, return a numeric vector with the number of 44 motifs: M111, M112, M113, M211, M212, M213, M311, M312, M121_1, M122_1, M122_2, M122_3, M123_1, M123_2, M123_3, M123_4, M123_5, M221_1, M221_2, M221_3, M222_1, M222_2, M222_3, M222_4, M222_5, M222_6, M222_7, M222_8, M321_1, M321_2, M321_3, M321_4, M321_5, M131, M132-1, M132-2, M132-3, M132-4, M132-5, M231-1, M231-2, M231-3, M231-4, M231-5.
#'
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
#' set.seed(12)
#' d <- build_net(11,22,21,0.2)
#' m <- motif_count(d)
#' m
#'
#' set.seed(12)
#' d <- build_net(11,22,21,0.2,asmatrices=TRUE)
#'
#' MAT<-d
#' M <- motif_count(MAT[[3]],MAT[[4]])
#' M
#'
#' R<-rownames(MAT[[4]])[12]
#' MR<-MAT[[4]][12,]
#' MAT[[4]]<-MAT[[4]][-12,]
#' MAT[[4]]<-rbind(MAT[[4]],MR)
#' rownames(MAT[[4]])[22]<-R
#'
#' motif_count(MAT[[3]],MAT[[4]])
#'
motif_count<-function(network.or.subnet_mat1,subnet_mat2=NULL,subnet_motif=FALSE){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network<-adject_net(network.or.subnet_mat1)
      PHP<-as.matrix(network[])
      dimnames(PHP)<-NULL
      PH<-PHP[(V(network)$level)==0,(V(network)$level)==1]
      HP<-PHP[(V(network)$level)==1,(V(network)$level)==2]
   }
   else if(inherits(network.or.subnet_mat1,c("matrix","data.frame"))==T && inherits(subnet_mat2,c("matrix","data.frame"))==T){
      mat1<-network.or.subnet_mat1
      mat2<-subnet_mat2
      if(ncol(mat1)!=nrow(mat2))
         stop("Error: please check whether the column of network.or.subnet_mat1 is corresponding to the row of subnet_mat2!!!")
      if(is.null(colnames(mat1)) || is.null(rownames(mat2))){
         colnames(mat1)<-paste0("mid_pse",seq=1:ncol(mat1))
         rownames(mat2)<-paste0("mid_pse",seq=1:ncol(mat1))
      }
      if(sum(!(colnames(mat1)%in%(rownames(mat2))),na.rm = TRUE)!=0 )
         stop("Error: please check whether the column name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
      if(sum(is.na(colnames(mat1)))!=0 || sum(is.na(rownames(mat2)))!=0)
         stop("Error: There is NA in the column name of network.or.subnet_mat1 or the row name of subnet_mat2!!!")
      mat2<-mat2[colnames(mat1),]
      dimnames(mat1)<-NULL
      dimnames(mat2)<-NULL
      PH<-mat1
      HP<-mat2
   }
   else
      stop("Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   L1<-ncol(PH)
   if(L1 < 4L)
      stop("Error: please input a large 'number of interconnecting species >=4' network data!!!")
   PH_add<-t(PH)%*%PH
   HP_add<-HP%*%t(HP)
   HP_ratio<-HP%*%(1-t(HP))
   PH_ratio<-t(PH)%*%(1-PH)
   Two<-function(s){
      sum(s)*(sum(s)-1)/2
   }
   Two_3<-function(s){
      sum(s)*(sum(s)-1)*(sum(s)-2)/6
   }
   PP<-PH%*%HP
   PHP_add<-PH_add*HP_add
   ############################################################
   M111<-sum(PP)##
   ############################################################
   M112=sum(PH%*%apply(HP,1,Two))
   ############################################################
   M211<-sum(apply((PH),2,Two)%*%HP)
   ############################################################
   M212<-sum(apply(PH,2,Two)%*%apply(HP,1,Two))
   ############################################################
   M121_1<-sum(PP*(PP-1)/2)##
   ############################################################
   M122_1<-sum(PH_add*HP_ratio*t(HP_ratio))/2
   ############################################################
   M122_2<-sum(PHP_add*HP_ratio)##
   ############################################################
   M<-PHP_add*(HP_add-1)##
   M122_3<-sum(M[lower.tri(M)])/2
   ############################################################
   M221_1<-sum(HP_add*PH_ratio*t(PH_ratio))/2
   ############################################################
   M221_2<-sum(PHP_add*(PH_ratio))##
   ############################################################
   M<-PHP_add*(PH_add-1)##
   M221_3<-sum(M[lower.tri(M)])/2
   ############################################################
   M222_1<-sum(PH_add*PH_ratio*HP_ratio*t(HP_ratio))
   ############################################################
   M222_2<-sum(PH_add*(PH_add-1)*HP_ratio*t(HP_ratio))/4
   ############################################################
   M222_3<-sum(PH_ratio*t(PH_ratio)*HP_add*HP_ratio)
   ############################################################
   M222_4<-sum(PHP_add*PH_ratio*HP_ratio)##
   ############################################################
   M222_5<-sum((PHP_add*(PH_add-1)/2)*HP_ratio)##
   ############################################################
   M222_6<-sum(PH_ratio*t(PH_ratio)*HP_add*(HP_add-1))/4
   ############################################################
   M222_7<-sum(PHP_add*PH_ratio*(HP_add-1)/2)##
   ############################################################
   M<-PHP_add*(PH_add-1)*(HP_add-1)/4##
   M222_8<-sum(M[lower.tri(M)])
   ###########################################################
   M222_9<-sum(PHP_add*PH_ratio*t(HP_ratio))##
   ###########################################################
   b<-apply(HP,1,Two_3)
   M113<-sum(PH%*%apply(HP,1,Two_3))
   ###########################################################
   M311<-sum(t(HP)%*%apply(t(PH),1,Two_3))
   ###########################################################
   M213<-sum(apply(PH,2,Two)%*%apply(HP, 1,Two_3))
   ###########################################################
   M312<-sum(apply(PH,2,Two_3)%*%apply(HP,1, Two))
   ###########################################################
   M123_1<-sum(PH_add*(HP_ratio*(HP_ratio-1)/2)*t(HP_ratio))
   ###########################################################
   M123_2<-sum(PHP_add*(HP_ratio*(HP_ratio-1)/2))##
   ###########################################################
   M123_3<-sum(PHP_add*HP_ratio*t(HP_ratio))/2##
   ###########################################################
   M123_4<-sum(PHP_add*((HP_add-1)/2)*HP_ratio)##
   ###########################################################
   M<-PHP_add*((HP_add-1)*(HP_add-2)/6)
   M123_5<-sum(M[lower.tri(M)])
   ###########################################################
   M321_1<-sum((PH_ratio*(PH_ratio-1)/2)*t(PH_ratio)*HP_add)
   ###########################################################
   M321_2<-sum((PH_ratio*(PH_ratio-1)/2)*PHP_add)##
   ###########################################################
   M321_3<-sum(PHP_add*PH_ratio*t(PH_ratio))/2##
   ############(##############################################
   M321_4<-sum((PHP_add*(PH_add-1)/2)*PH_ratio)##
   ###########################################################
   M<-(PHP_add*(PH_add-1)*(PH_add-2)/6)##
   M321_5<-sum(M[lower.tri(M)])
   ###########################################################
   M131<-sum(PP*(PP-1)*(PP-2)/6)##
   ###########################################################
   M1321=M1322=M1323=M1324=M1325=M2311=M2312=M2313=M2314=M2315<-0
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
   motif<-c(M111,M112,M113,M211,M212,M213,M311,M312,M121_1,M122_1,M122_2,M122_3,M123_1,M123_2,M123_3,M123_4,M123_5,M221_1,M221_2,M221_3,M222_1,M222_2,M222_3,M222_4,M222_5,M222_6,M222_7,M222_8,M321_1,M321_2,M321_3,M321_4,M321_5,M131,M1321,M1322,M1323,M1324,M1325,M2311,M2312,M2313,M2314,M2315)
   ###########################################################
   # if(subnet_motif){
   #    subnet1_motif<-bmotif::mcount(PH,six_node = T, normalisation = TRUE, mean_weight = F, standard_dev = F)
   #    subnet2_motif<-bmotif::mcount(HP,six_node = T, normalisation = TRUE, mean_weight = F, standard_dev = F)
   #    subnet1_node_position<-bmotif::node_positions(PH,six_node=TRUE,weights_method="none")
   #    subnet2_node_position<-bmotif::node_positions(HP,six_node=TRUE,weights_method="none")
   #    subnet1_link_position<-bmotif::link_positions(PH,six_node=TRUE,weights =FALSE)
   #    subnet2_link_position<-bmotif::link_positions(HP,six_node=TRUE,weights =FALSE)
   #    motif_list<-list()
   #    motif_list$multilayernet_motif<-motif
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
