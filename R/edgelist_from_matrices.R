edgelist_from_matrices<-function(mat1,mat2,isDirected1=T,isDirected2=T){
   if(ncol(mat1)!=nrow(mat2))
      stop("Error: please check whether the column of mat1 is corresponding to the row of mat2!!!")
   if(is.null(colnames(mat1)) || is.null(rownames(mat2))){
      colnames(mat1)<-paste0("mid_pse",seq=1:ncol(mat1))
      rownames(mat2)<-paste0("mid_pse",seq=1:ncol(mat1))
   }
   if(sum(!(colnames(mat1)%in%(rownames(mat2))),na.rm = TRUE)!=0 )
      stop("Error: please check whether the column name of mat1 is corresponding to the row name of mat2!!!")
   if(sum(is.na(colnames(mat1)))!=0 || sum(is.na(rownames(mat2)))!=0)
      stop("Error: There is NA in the column name of mat1 or the row name of mat2!!!")
   mat2<-mat2[colnames(mat1),]
   if(is.null(rownames(mat1)))
      rownames(mat1)<-paste0("up_pse",seq=1:nrow(mat1))
   if(is.null(colnames(mat2)))
      colnames(mat2)<-paste0("down_pse",seq=1:ncol(mat2))
   spe<-unique(c(rownames(mat1),colnames(mat1),rownames(mat2),colnames(mat2)))
   me_interlayer<-function(mat,U,V,isDirected=FALSE){
      spe_v<-NULL
      mat[mat[]>0]<-1
      mat[is.na(mat)]<-0
      for (i in 1:ncol(mat)) {
         spe_v<-c(spe_v,mat[,i])
      }
      spe_mat<-data.frame(spe1=rep(rownames(mat),ncol(mat)),layer=rep(U,ncol(mat)*nrow(mat)),spe2=(rep(colnames(mat),each=nrow(mat))),layer1=rep(V,nrow(mat)*ncol(mat)),value=spe_v)
      colnames(spe_mat)<-NULL
      spe_mat<-as.matrix(spe_mat)
      spe_mat<-spe_mat[spe_mat[,5]==1,]
      if(isDirected==FALSE){
         spe_m<-spe_mat
         spe_m[,1:2]<-spe_mat[,3:4]
         spe_m[,3:4]<-spe_mat[,1:2]
         spe_m<-as.matrix(spe_m)
         spe_mat<-as.matrix(spe_mat)
         spe_mat<-rbind(spe_mat,spe_m)
      }
      rownames(spe_mat)<-NULL
      return(spe_mat)
   }
   matt<-rbind(me_interlayer(mat1,1,2,isDirected1),me_interlayer(mat2,2,3,isDirected2))
   for (i in 1:sum(!(is.na(spe)))) {
      matt[matt[,1]%in%spe[i],1]<-i
      matt[matt[,3]%in%spe[i],3]<-i
   }
   matt<-as.data.frame(matt)
   matt<-apply(matt, 2,function(x){as.numeric(x)})
   degree<-c(rowSums(mat1),colSums(mat1)+rowSums(mat2),colSums(mat2))
   return(list(matt,spe=spe,degree=degree))
}
