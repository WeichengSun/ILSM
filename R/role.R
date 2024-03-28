#' Analyzising role of interconnecting node in motifs
#'
#' Counting the number of 65 roles about interconnecting species in multilayer network motifs.
#'
#' @param network.or.subnet_mat1 Either a multilayer(tripartite) network of 'igraph' class, or a numeric matrix(or data.frame) representing interactions between two groups of species. The network contains
#'  interlayer links and without intralayer links. Each row and column of matrix represents the species in the first and second layers of the tripartite network respectively.
#'  Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise.
#'
#' @param subnet_mat2 A numeric matrix(or data.frame) representing interactions between two groups of species.Each row and column of matrix represents the species in the second and third layers of the
#'  tripartite network respectively. Elements of matrix are non-zero numebers if the interlayer species are connected, and 0 otherwise. If \code{network.or.subnet_mat1} is "igraph", \code{subnet_mat2} defaults to NULL.
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
#' @import igraph
#' @export
#'
#' @return
#' Returns a matrix of 65 columns representing the roles of interconnecting species in the motifs. Columns names are Role1, Role2, Role3 ... Role65.
#'
#' Each row of matrix corresponds to a interconnecting species in the second layer of network. If a interconnecting species is linked to both the second and third level species, the elements in this row are not all zero, otherwise the elements are all zero.
#'
#' @references
#' Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V., Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif analyses of bipartite networks. Methods in Ecology and Evolution, 10(5), 695-701.
#'
#' @examples
#'
#' set.seed(12)
#' d <- build_net(11,22,21,0.2)
#' mr <- Midlayer_role(d)
#' mr
#'
#' set.seed(12)
#' MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
#'
#' M <- Midlayer_role(MAT[[3]],MAT[[4]])
#' M
#'
#' R<-rownames(MAT[[4]])[12]
#' MR<-MAT[[4]][12,]
#' MAT[[4]]<-MAT[[4]][-12,]
#' MAT[[4]]<-rbind(MAT[[4]],MR)
#' rownames(MAT[[4]])[22]<-R
#'
#' Midlayer_role(MAT[[3]],MAT[[4]])
#'
Midlayer_role<-function(network.or.subnet_mat1,subnet_mat2=NULL){
   if(inherits(network.or.subnet_mat1,"igraph")==T){
      network<-adject_net(network.or.subnet_mat1)
      PHP<-as.matrix(network[])
      PH<-PHP[(V(network)$level)==0,(V(network)$level)==1]
      HP<-PHP[(V(network)$level)==1,(V(network)$level)==2]
      spe<-V(network.or.subnet_mat1)$name[V(network.or.subnet_mat1)$level==1]
      role<-matrix(0,length(spe),65)
      rownames(role)<-spe
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
      PH<-mat1
      HP<-mat2
      logi<-(apply(PH,2,sum)*apply(HP,1,sum))!=0
      PH<-PH[,logi]
      HP<-HP[logi,]
      spe<-colnames(mat1)[colnames(mat1)==rownames(mat2)]
      role<-matrix(0,length(spe),65)
      rownames(role)<-spe
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
   HH_role<-NULL

   M111<-(apply(PH,2,sum)*apply(HP,1,sum))
   HH_role<-cbind(HH_role,M111)

   M112<-(apply(PH,2,sum)*apply(HP,1,function(i){sum(i)*(sum(i)-1)/2}))
   HH_role<-cbind(HH_role,M112)

   M113<-(apply(PH,2,sum)*apply(HP,1,function(i){sum(i)*(sum(i)-1)*(sum(i)-2)/6}))
   HH_role<-cbind(HH_role,M113)

   M211<-apply(PH,2,function(i){sum(i)*(sum(i)-1)/2})*apply(HP,1,sum)
   HH_role<-cbind(HH_role,M211)

   M212<-(apply(PH,2,function(i){sum(i)*(sum(i)-1)/2})*apply(HP,1,function(i){sum(i)*(sum(i)-1)/2}))
   HH_role<-cbind(HH_role,M212)

   M213<-apply(PH,2,function(i){sum(i)*(sum(i)-1)/2})*apply(HP,1,function(i){sum(i)*(sum(i)-1)*(sum(i)-2)/6})
   HH_role<-cbind(HH_role,M213)

   M311<-apply(PH,2,function(i){sum(i)*(sum(i)-1)*(sum(i)-2)/6})*apply(HP,1,sum)
   HH_role<-cbind(HH_role,M311)

   M312<-apply(PH,2,function(i){sum(i)*(sum(i)-1)*(sum(i)-2)/6})*apply(HP,1,function(i){sum(i)*(sum(i)-1)/2})
   HH_role<-cbind(HH_role,M312)

   M121_1<-(PH_add*HP_add)*((PH_add*HP_add)%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M121_1)+colSums(M121_1))

   M122_1<-(PH_add*HP_ratio*t(HP_ratio))*((PH_add*HP_ratio*t(HP_ratio))%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M122_1)+colSums(M122_1))

   M122_2<-PH_add*HP_add*HP_ratio
   HH_role<-cbind(HH_role,colSums(M122_2),rowSums(M122_2))

   M122_3<-(PH_add*HP_add*(HP_add-1)/2)*((PH_add*HP_add*(HP_add-1)/2)%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M122_3)+colSums(M122_3))

   M123_1<-PH_add*(HP_ratio*(HP_ratio-1)/2)*t(HP_ratio)
   HH_role<-cbind(HH_role,rowSums(M123_1),colSums(M123_1))

   M123_2<-PH_add*HP_add*(HP_ratio*(HP_ratio-1)/2)
   HH_role<-cbind(HH_role,colSums(M123_2),rowSums(M123_2))

   M123_3<-PH_add*HP_add*HP_ratio*t(HP_ratio)/2
   HH_role<-cbind(HH_role,rowSums(M123_3)+colSums(M123_3))

   M123_4<-PH_add*(HP_add*(HP_add-1)/2)*HP_ratio
   HH_role<-cbind(HH_role,rowSums(M123_4),colSums(M123_4))

   M123_5<-(PH_add*(HP_add*(HP_add-1)*(HP_add-2)/6))*((PH_add*(HP_add*(HP_add-1)*(HP_add-2)/6))%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M123_5)+colSums(M123_5))

   M221_1<-(HP_add*PH_ratio*t(PH_ratio))* ((HP_add*PH_ratio*t(PH_ratio))%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M221_1)+colSums(M221_1))

   M221_2<-PH_add*(PH_ratio)*HP_add
   HH_role<-cbind(HH_role,rowSums(M221_2),colSums(M221_2))

   M221_3<-(PH_add*(PH_add-1)*HP_add/2)*((PH_add*(PH_add-1)*HP_add/2)%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M221_3)+colSums(M221_3))

   M222_1<-PH_add*PH_ratio*HP_ratio*t(HP_ratio)
   HH_role<-cbind(HH_role,rowSums(M222_1),colSums(M222_1))

   M222_2<-PH_add*(PH_add-1)*HP_ratio*t(HP_ratio)/4
   HH_role<-cbind(HH_role,rowSums(M222_2)+colSums(M222_2))

   M222_3<-PH_ratio*t(PH_ratio)*HP_add*HP_ratio
   HH_role<-cbind(HH_role,rowSums(M222_3),colSums(M222_3))

   M222_4<-PH_add*HP_add*PH_ratio*HP_ratio
   HH_role<-cbind(HH_role,rowSums(M222_4),colSums(M222_4))

   M222_5<-(PH_add*(PH_add-1)/2)*HP_add*HP_ratio
   HH_role<-cbind(HH_role,colSums(M222_5),rowSums(M222_5))

   M222_6<-PH_ratio*t(PH_ratio)*HP_add*(HP_add-1)/4
   HH_role<-cbind(HH_role,rowSums(M222_6)+colSums(M222_6))

   M222_7<-PH_add*PH_ratio*HP_add*(HP_add-1)/2
   HH_role<-cbind(HH_role,rowSums(M222_7),colSums(M222_7))

   M222_8<-(PH_add*(PH_add-1)*HP_add*(HP_add-1)/4)*((PH_add*(PH_add-1)*HP_add*(HP_add-1)/4)%>%upper.tri())
   HH_role<-cbind(HH_role,rowSums(M222_8)+colSums(M222_8))

   M321_1<-(PH_ratio*(PH_ratio-1)/2)*t(PH_ratio)*HP_add
   HH_role<-cbind(HH_role,rowSums(M321_1),colSums(M321_1))

   M321_2<-(PH_ratio*(PH_ratio-1)/2)*PH_add*HP_add
   HH_role<-cbind(HH_role,colSums(M321_2),rowSums(M321_2))

   M321_3<-(PH_add*PH_ratio*t(PH_ratio)*HP_add/2)
   HH_role<-cbind(HH_role,colSums(M321_3)+rowSums(M321_3))

   M321_4<-(PH_add*(PH_add-1)/2)*PH_ratio*HP_add
   HH_role<-cbind(HH_role,colSums(M321_4),rowSums(M321_4))

   M321_5<-(HP_add*PH_add*(PH_add-1)*(PH_add-2)/6)*((HP_add*PH_add*(PH_add-1)*(PH_add-2)/6)%>%upper.tri())
   HH_role<-cbind(HH_role,colSums(M321_5)+rowSums(M321_5))

####################
   M131=M1321=M1321_1=M1322=M1322_1=M1323=M1323_1=M1324=M1324_1=M1325=M2311=M2311_1=M2312=M2312_1=M2313=M2313_1=M2314=M2314_1=M2315<-rep(0,L1)
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
      }
   }
   HH_role<-cbind(HH_role,M131/6,M1321/2,M1321_1/2,M1322_1/2,M1322/2,M1323/2,M1323_1/2,M1324/2,M1324_1/2,M1325/6,M2311/2,M2311_1/2,M2312/2,M2312_1/2,M2313/2,M2313_1/2,M2314/2,M2314_1/2,M2315/6)
   colnames(HH_role)<-paste0("role",c(1:65))
   role[rownames(HH_role),]<-HH_role
   colnames(role)<-paste0("role",c(1:65))
   return(role)
}
