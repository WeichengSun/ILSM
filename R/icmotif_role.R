#' Analyzing role of interconnecting node in motifs
#'
#' Counting the number of 70 roles about interconnecting species in tripartite network motifs.
#'
#' @param network.or.subnet_mat1 Either a tripartite network of 'igraph' class which contains three groups of species and interactions within subnetwork without interactions between each group of species, or a numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and first groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.
#'  Each row and column of matrix represents single species in the second and third groups of the tripartite network respectively.
#'  Elements of matrix are non-zero numbers if the two groups of species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
#'
#' @param ic_number The max number of interconnection code in each motif. Defaults to 4 because of containing all motifs, also can input 2 to calculate the top 36 for rapid computing speed.
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
#' @import igraph
#' @export
#'
#' @return
#' If \code{ic_number} = 4,
#' return a matrix of 70 columns representing the roles of interconnecting
#' species in the motifs. Columns names are Role1, Role2, Role3 ... Role70.
#'
#' Wherever, the \code{ic_number} = 2,
#' return a matrix of 50 columns representing the roles of interconnecting
#' species in the top 30 motifs. Columns names are Role1, Role2,
#' Role3 ... Role50.
#'
#' Each row of matrix corresponds to a interconnecting species in the second
#' group of network. If a interconnecting species is linked to both the first
#' and third group species, the elements in this row are not all zero, otherwise
#'  the elements are all zero.
#'
#' @srrstats {G1.1} The algorithm is the first implementation of a novel
#' lgorithm.
#'
#' @references
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#' @examples
#'
#' ## generate a random tripartite network
#' set.seed(12)
#' Net <- build_net(11,15,16,0.2)
#'
#' data(PPH_Coltparkmeadow)
#' Net <- PPH_Coltparkmeadow
#' icmotif_role(Net)
#' icmotif_role(Net, ic_number=2)
#'
#' set.seed(12)
#' MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
#'
#' icmotif_role(MAT[[3]],MAT[[4]])
#'
#' md1 <- matrix(sample(c(0,1),88,replace=TRUE),8,11)
#' md2 <- matrix(sample(c(0,1),120,replace=TRUE),10,12)
#' icmotif_role(md1,md2)
#'
#' R <- rownames(MAT[[4]])[12]
#' MR <- MAT[[4]][12,]
#' MAT[[4]] <- MAT[[4]][-12,]
#' MAT[[4]] <- rbind(MAT[[4]],MR)
#' rownames(MAT[[4]])[22] <- R
#'
#' icmotif_role(MAT[[3]],MAT[[4]])
#'
#'

icmotif_role<-function(network.or.subnet_mat1, subnet_mat2=NULL, ic_number=4){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network<-adjust_net(network.or.subnet_mat1)
      PHP<-as.matrix(network[])
      PHP[PHP>0]<-1
      PH<-PHP[(V(network)$level)==0,(V(network)$level)==1]
      HP<-PHP[(V(network)$level)==1,(V(network)$level)==2]
      spe<-V(network.or.subnet_mat1)$name[V(network.or.subnet_mat1)$level==1]
      role<-matrix(0,length(spe),70)
      rownames(role)<-spe
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

      PH<-t(mat1)
      HP<-mat2
      logi<-(apply(PH,2,sum)*apply(HP,1,sum))!=0
      PH<-PH[,logi]
      HP<-HP[logi,]
      spe<-rownames(mat1)[rownames(mat1)==rownames(mat2)]
      role<-matrix(0,length(spe),70)
      rownames(role)<-spe
   }
   else
      stop("Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   L1<-ncol(PH)
   if(L1 < 4L)
      stop("Error: please input a large 'number of interconnecting species >=4' network data!!!")
   # PH_add<-t(PH)%*%PH
   # HP_add<-HP%*%t(HP)
   # HP_ratio<-HP%*%(1-t(HP))
   # PH_ratio<-t(PH)%*%(1-PH)
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
   HH_role<-NULL

   M111 <- Ob*Rb
   HH_role <- cbind(HH_role,M111)

   M112 <- Ob*Two(Rb)
   HH_role <- cbind(HH_role,M112)

   M113 <- Ob*Three(Rb)
   HH_role <- cbind(HH_role,M113)

   M114<-Ob*Four(Rb)
   HH_role <- cbind(HH_role,M114)

   M211 <- Two(Ob)*Rb
   HH_role <- cbind(HH_role,M211)

   M212 <- Two(Ob)*Two(Rb)
   HH_role <- cbind(HH_role,M212)

   M213 <- Two(Ob)*Three(Rb)
   HH_role <- cbind(HH_role,M213)

   M311 <- Three(Ob)*Rb
   HH_role <- cbind(HH_role,M311)

   M312 <- Three(Ob)*Two(Rb)
   HH_role <- cbind(HH_role,M312)

   M411 <- Four(Ob)*Rb
   HH_role <- cbind(HH_role,M411)

   # M121#
   HH_role <- cbind(HH_role,rowSums(UVbb)-M111)

   M122_1 <- Ubb*Ybb*YT
   HH_role <- cbind(HH_role,rowSums(M122_1))

   M122_2 <- UVbb*Ybb
   HH_role <- cbind(HH_role,colSums(M122_2),rowSums(M122_2))

   M122_3 <- Ubb*Two(Vbb)
   HH_role <- cbind(HH_role,rowSums(M122_3)-M112)

   M123_1 <- Ubb*Two(Ybb)*YT
   HH_role <- cbind(HH_role,rowSums(M123_1),colSums(M123_1))

   M123_2 <- UVbb*Two(YT)
   HH_role <- cbind(HH_role,rowSums(M123_2),colSums(M123_2))

   M123_3 <- UVbb*Ybb*YT
   HH_role <- cbind(HH_role,rowSums(M123_3))

   M123_4 <- Ubb*Two(Vbb)*Ybb
   HH_role <- cbind(HH_role,rowSums(M123_4),colSums(M123_4))

   M123_5 <- Ubb*Three(Vbb)
   HH_role <- cbind(HH_role,rowSums(M123_5)-M113)

   M221_1 <- Xbb*XT*Vbb
   HH_role <- cbind(HH_role,rowSums(M221_1))

   M221_2 <- UVbb*Xbb
   HH_role <- cbind(HH_role,rowSums(M221_2),colSums(M221_2))

   M221_3 <- Two(Ubb)*Vbb
   HH_role <- cbind(HH_role,rowSums(M221_3)-M211)

   M222_1 <- Ubb*Xbb*Ybb*YT
   HH_role <- cbind(HH_role,rowSums(M222_1),colSums(M222_1))

   M222_2 <- Two(Ubb)*Ybb*YT
   HH_role <- cbind(HH_role,rowSums(M222_2))

   M222_3 <- Xbb*XT*Ybb*Vbb
   HH_role <- cbind(HH_role,rowSums(M222_3),colSums(M222_3))

   M222_4 <- UVbb*Xbb*Ybb
   HH_role <- cbind(HH_role,rowSums(M222_4),colSums(M222_4))

   M222_5 <- Two(Ubb)*Vbb*YT
   HH_role <- cbind(HH_role,rowSums(M222_5),colSums(M222_5))

   M222_6 <- Xbb*XT*Two(Vbb)
   HH_role <- cbind(HH_role,rowSums(M222_6))

   M222_7 <- Ubb*Xbb*Two(Vbb)
   HH_role <- cbind(HH_role,rowSums(M222_7),colSums(M222_7))

   M222_8 <- Two(Ubb)*Two(Vbb)
   HH_role <- cbind(HH_role,rowSums(M222_8)-M212)

   M222_9 <- UVbb*Xbb*YT
   HH_role <- cbind(HH_role,rowSums(M222_9),colSums(M222_9))

   M321_1 <- Two(Xbb)*XT*Vbb
   HH_role <- cbind(HH_role,rowSums(M321_1),colSums(M321_1))

   M321_2 <- UVbb*Two(XT)
   HH_role <- cbind(HH_role,rowSums(M321_2),colSums(M321_2))

   M321_3 <- UVbb*Xbb*XT
   HH_role <- cbind(HH_role,rowSums(M321_3))

   M321_4 <- Two(Ubb)*XT*Vbb
   HH_role <- cbind(HH_role,rowSums(M321_4),colSums(M321_4))

   M321_5 <- Three(Ubb)*Vbb
   HH_role <- cbind(HH_role,rowSums(M321_5)-M311)

   if(ic_number==2){
      role<-role[,1:50]
      role[rownames(HH_role),]<-HH_role
      colnames(role)<-paste0("role",c(1:50))
      return(role)
   }

####################
   M131<-M1321<-M1321_1<-M1322<-M1322_1<-M1323<-M1323_1<-M1324<-rep(0,L1)
   M1324_1<-M1325<-M2311<-M2311_1<-M2312<-M2312_1<-M2313<-M2313_1<-rep(0,L1)
   M2314=M2314_1=M2315=M141<-rep(0,L1)
   if(L1>4){
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
               role_value<-PH_F*HP_F
               if(sum(role_value)!=0){
                  M131[c(i,j)]<-M131[c(i,j)]+sum(role_value)
                  M131[-c(i,j)]<-M131[-c(i,j)]+role_value
               }
               role_value1<-PH_F*HP_F1*HP_F2
               if(sum(role_value1)!=0){
                  M1321[c(i,j)]<-M1321[c(i,j)]+sum(role_value1)
                  M1321_1[-c(i,j)]<-M1321_1[-c(i,j)]+role_value1
               }
               role_value2<-PH_F*HP_F*HP_F2
               if(sum(role_value2)!=0){
                  M1322[c(i,j)]<-M1322[c(i,j)]+sum(role_value2)
                  M1322_1[-c(i,j)]<-M1322_1[-c(i,j)]+role_value2
               }
               role_value3<-PH_F*HP_F1*HP_F3
               if(sum(role_value3)!=0){
                  M1323[i]<-M1323[i]+sum(role_value3)
                  M1323[-c(i,j)]<-M1323[-c(i,j)]+role_value3
                  M1323_1[j]<-M1323_1[j]+sum(role_value3)
               }
               role_value4<-PH_F*HP_F*HP_F1
               if(sum(role_value4)!=0){
                  M1324[c(i,j)]<-M1324[c(i,j)]+sum(role_value4)
                  M1324_1[-c(i,j)]<-M1324_1[-c(i,j)]+role_value4
               }
               role_value5<-PH_F*HP_F*(HP_F-1)/2
               if(sum(role_value5)!=0){
                  M1325[c(i,j)]<-M1325[c(i,j)]+sum(role_value5)
                  M1325[-c(i,j)]<-M1325[-c(i,j)]+role_value5
               }
               role_value6<-PH_F1*PH_F2*HP_F
               if(sum(role_value6)!=0){
                  M2311[c(i,j)]<-M2311[c(i,j)]+sum(role_value6)
                  M2311_1[-c(i,j)]<-M2311_1[-c(i,j)]+role_value6
               }
               role_value7<-PH_F*PH_F2*HP_F
               if(sum(role_value7)!=0){
                  M2312_1[c(i,j)]<-M2312_1[c(i,j)]+sum(role_value7)
                  M2312[-c(i,j)]<-M2312[-c(i,j)]+(role_value7)
               }
               role_value8<-PH_F1*PH_F3*HP_F
               if(sum(role_value8)!=0){
                  M2313[i]<-M2313[i]+sum(role_value8)
                  M2313[-c(i,j)]<-M2313[-c(i,j)]+role_value8
                  M2313_1[j]<-M2313_1[j]+sum(role_value8)
               }
               role_value9<-PH_F*PH_F3*HP_F
               if(sum(role_value9)!=0){
                  M2314[i]<-M2314[i]+sum(role_value9)
                  M2314_1[j]<-M2314_1[j]+sum(role_value9)
                  M2314_1[-c(i,j)]<-M2314_1[-c(i,j)]+role_value9
               }
               role_value10<-PH_F*(PH_F-1)*HP_F/2
               if(sum(role_value10)!=0){
                  M2315[c(i,j)]<-M2315[c(i,j)]+sum(role_value10)
                  M2315[-c(i,j)]<-M2315[-c(i,j)]+role_value10
               }
            }
            if(i<j){
               for(k in (j):L1){
                  if(j<k){
                     role_value11<-((PH[,i]*PH[,j]*PH[,k])%*%PH[,-c(i,j,k)])*((HP[i,]*HP[j,]*HP[k,])%*%t(HP[-c(i,j,k),]))
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
               PH_F<-colSums(PH[,i]*PH[,j]*(PH[,-c(i,j)]))
               PH_F1<-colSums(PH[,i]*PH[,j]*(1-PH[,-c(i,j)]))
               PH_F2<-colSums((1-PH[,i])*(1-PH[,j])*(PH[,-c(i,j)]))
               PH_F3<-colSums((1-PH[,i])*(PH[,j])*(PH[,-c(i,j)]))
               HP_F<-colSums(HP[i,]*HP[j,]*t(HP[-c(i,j),]))
               HP_F1<-colSums(HP[i,]*HP[j,]*t(1-HP[-c(i,j),]))
               HP_F2<-colSums((1-HP[i,])*(1-HP[j,])*t(HP[-c(i,j),]))
               HP_F3<-colSums((1-HP[i,])*(HP[j,])*t(HP[-c(i,j),]))
               ########
               role_value<-PH_F*HP_F
               if(sum(role_value)!=0){
                  M131[c(i,j)]<-M131[c(i,j)]+sum(role_value)
                  M131[-c(i,j)]<-M131[-c(i,j)]+role_value
               }
               role_value1<-PH_F*HP_F1*HP_F2
               if(sum(role_value1)!=0){
                  M1321[c(i,j)]<-M1321[c(i,j)]+sum(role_value1)
                  M1321_1[-c(i,j)]<-M1321_1[-c(i,j)]+role_value1
               }
               role_value2<-PH_F*HP_F*HP_F2
               if(sum(role_value2)!=0){
                  M1322[c(i,j)]<-M1322[c(i,j)]+sum(role_value2)
                  M1322_1[-c(i,j)]<-M1322_1[-c(i,j)]+role_value2
               }
               role_value3<-PH_F*HP_F1*HP_F3
               if(sum(role_value3)!=0){
                  M1323[i]<-M1323[i]+sum(role_value3)
                  M1323[-c(i,j)]<-M1323[-c(i,j)]+role_value3
                  M1323_1[j]<-M1323_1[j]+sum(role_value3)
               }
               role_value4<-PH_F*HP_F*HP_F1
               if(sum(role_value4)!=0){
                  M1324[c(i,j)]<-M1324[c(i,j)]+sum(role_value4)
                  M1324_1[-c(i,j)]<-M1324_1[-c(i,j)]+role_value4
               }
               role_value5<-PH_F*HP_F*(HP_F-1)/2
               if(sum(role_value5)!=0){
                  M1325[c(i,j)]<-M1325[c(i,j)]+sum(role_value5)
                  M1325[-c(i,j)]<-M1325[-c(i,j)]+role_value5
               }
               role_value6<-PH_F1*PH_F2*HP_F
               if(sum(role_value6)!=0){
                  M2311[c(i,j)]<-M2311[c(i,j)]+sum(role_value6)
                  M2311_1[-c(i,j)]<-M2311_1[-c(i,j)]+role_value6
               }
               role_value7<-PH_F*PH_F2*HP_F
               if(sum(role_value7)!=0){
                  M2312_1[c(i,j)]<-M2312_1[c(i,j)]+sum(role_value7)
                  M2312[-c(i,j)]<-M2312[-c(i,j)]+(role_value7)
               }
               role_value8<-PH_F1*PH_F3*HP_F
               if(sum(role_value8)!=0){
                  M2313[i]<-M2313[i]+sum(role_value8)
                  M2313[-c(i,j)]<-M2313[-c(i,j)]+role_value8
                  M2313_1[j]<-M2313_1[j]+sum(role_value8)
               }
               role_value9<-PH_F*PH_F3*HP_F
               if(sum(role_value9)!=0){
                  M2314[i]<-M2314[i]+sum(role_value9)
                  M2314_1[j]<-M2314_1[j]+sum(role_value9)
                  M2314_1[-c(i,j)]<-M2314_1[-c(i,j)]+role_value9
               }
               role_value10<-PH_F*(PH_F-1)*HP_F/2
               if(sum(role_value10)!=0){
                  M2315[c(i,j)]<-M2315[c(i,j)]+sum(role_value10)
                  M2315[-c(i,j)]<-M2315[-c(i,j)]+role_value10
               }
            }
            M141<-sum(PH[,1]*PH[,2]*PH[,3]*PH[,4])*sum(HP[1,]*HP[2,]*HP[3,]*HP[4,])
         }
      }
   }
   HH_role<-cbind(HH_role,M131/6,M1321/2,M1321_1/2,M1322_1/2,M1322/2,M1323/2,M1323_1/2,M1324/2,M1324_1/2,M1325/6,M2311/2,M2311_1/2,M2312/2,M2312_1/2,M2313/2,M2313_1/2,M2314/2,M2314_1/2,M2315/6,M141)
   colnames(HH_role)<-paste0("role",c(1:70))
   role[rownames(HH_role),]<-HH_role
   colnames(role)<-paste0("role",c(1:70))
   return(role)
}
