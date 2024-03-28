#' @importFrom muxViz BuildSupraAdjacencyMatrixFromExtendedEdgelist
#' @importFrom muxViz GetMultiPageRankCentrality
#' @importFrom muxViz GetMultiDegreeSum
#' @importFrom muxViz GetMultiHubCentrality
#' @importFrom muxViz GetMultiAuthCentrality
#' @importFrom muxViz GetMultiKatzCentrality
#' @importFrom muxViz GetMultiEigenvectorCentrality
#' @importFrom muxViz GetMultiClosenessCentrality

Node_versatility1 <- function(mat1,mat2,mat3,type=c("degree","pagerank","multiplexity","hub","authority","katz","kcore","eigenvector","closeness","all")){
   data_trans<-function(mat1,mat2,mat3){
      spe<-unique(c(rownames(mat1),colnames(mat1),rownames(mat2),colnames(mat2),rownames(mat3),colnames(mat3)))

      me<-function(mat,U){
         spe_v<-NULL
         mat[mat[]>0]<-1
         mat[is.na(mat)]<-0
         for (i in 1:ncol(mat)) {
            spe_v<-c(spe_v,mat[,i])
         }
         spe_mat<-data.frame(spe1=rep(rownames(mat),ncol(mat)),layer=rep(U,ncol(mat)*nrow(mat)),spe2=(rep(colnames(mat),each=nrow(mat))),layer1=rep(U,nrow(mat)*ncol(mat)),value=spe_v)
         colnames(spe_mat)<-NULL

         spe_mat<-spe_mat[spe_mat[,5]==1,]

         spe_m<-spe_mat
         spe_m[,1]<-spe_mat[,3]
         spe_m[,3]<-spe_mat[,1]

         spe_m<-as.matrix(spe_m)
         spe_mat<-as.matrix(spe_mat)
         spe_mat<-rbind(spe_mat,spe_m)
         rownames(spe_mat)<-NULL
         return(spe_mat)
      }
      matt<-rbind(me(mat1,1),me(mat2,2),me(mat3,3))
      spe1<-unique(c(rownames(mat1),colnames(mat1)))
      spe2<-unique(c(rownames(mat2),colnames(mat2)))
      spe3<-unique(c(rownames(mat3),colnames(mat3)))
      spe12<-intersect(spe1,spe2)

      spe23<-intersect(spe2,spe3)

      spe12m<-data.frame(s1=rep(spe12,each=2),layer=rep(c(1,2),length(spe12)),s2=rep(spe12,each=2),layer1=rep(c(2,1),length(spe12)),value=1)%>%as.matrix();colnames(spe12m)<-NULL
      spe23m<-data.frame(s1=rep(spe23,each=2),layer=rep(c(2,3),length(spe23)),s2=rep(spe23,each=2),layer1=rep(c(3,2),length(spe23)),value=1)%>%as.matrix();colnames(spe23m)<-NULL
      matt<-rbind(matt,spe12m,spe23m)

      for (i in 1:sum(!(is.na(spe)))) {
         matt[matt[,1]%in%spe[i],1]<-i
         matt[matt[,3]%in%spe[i],3]<-i
      }
      matt<-as.data.frame(matt)
      matt[,1]<-as.numeric(matt[,1])
      matt[,2]<-as.numeric(matt[,2])
      matt[,3]<-as.numeric(matt[,3])
      matt[,4]<-as.numeric(matt[,4])
      matt[,5]<-as.numeric(matt[,5])
      # matt<-apply(matt, 2,function(x){as.matrix(x)})
      return(list(matt,spe=spe))

   }
   EDGES <- data_trans(mat1,mat2,mat3)[[1]]
   node <- data_trans(mat1,mat2,mat3)[[2]]
   NODE <- length(data_trans(mat1,mat2,mat3)[[2]])
   LAYER <- 3
   SupraAdjacencyMatrix <- BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges=EDGES, Layers=LAYER, Nodes=NODE, isDirected=FALSE)
   if(missing(type))
      type <- "all"
   if (type == "degree")
      return(data.frame(node=node, degree=GetMultiDegreeSum(SupraAdjacencyMatrix, LAYER, NODE, isDirected=F)))
   if (type == "pagerank")
      return(data.frame(node=node, PageRank_versatility=GetMultiPageRankCentrality(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "hub")
      return(data.frame(node=node, Hub_versatility=GetMultiHubCentrality(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "authority")
      return(data.frame(node=node, Authority_versatility=GetMultiAuthCentrality(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "katz")
      return(data.frame(node=node, Katz_versatility=GetMultiKatzCentrality(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "eigenvector")
      return(data.frame(node=node, Eigenvector_versatility=GetMultiEigenvectorCentrality(SupraAdjacencyMatrix, LAYER, NODE)))
   if (type == "closeness")
      return(data.frame(node=node, Closeness_versatility=GetMultiClosenessCentrality(SupraAdjacencyMatrix, LAYER, NODE)))

   if(type == "all"){
      degree <- GetMultiDegreeSum(SupraAdjacencyMatrix, LAYER, NODE, isDirected=F)
      PageRank_versatility <- GetMultiPageRankCentrality(SupraAdjacencyMatrix, LAYER, NODE)
      Hub_versatility <- GetMultiHubCentrality(SupraAdjacencyMatrix, LAYER, NODE)
      Authority_versatility <- GetMultiAuthCentrality(SupraAdjacencyMatrix, LAYER, NODE)
      Katz_versatility <- GetMultiKatzCentrality(SupraAdjacencyMatrix, LAYER, NODE)
      Eigenvector_versatility <- GetMultiEigenvectorCentrality(SupraAdjacencyMatrix, LAYER, NODE)
      Closeness_versatility <- GetMultiClosenessCentrality(SupraAdjacencyMatrix, LAYER, NODE)
      versatility <- data.frame(node=node, degree=degree, PageRank_versatility=PageRank_versatility, Hub_versatility=Hub_versatility, Authority_versatility=Authority_versatility,
                                Katz_versatility=Katz_versatility, Eigenvector_versatility=Eigenvector_versatility, Closeness_versatility=Closeness_versatility)
      return(versatility)
   }
   else
      stop("Error: type is not a valid input!")
}
