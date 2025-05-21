SauveR<-function(N,m,type=c("row","col")){
   M<-emptyR(m)
   web <- as.matrix(M[[1]])
   if(type=="row"){
      sauve.fastr <- function(web,m,M) {
         web[sample(sample(nrow(web))),]<-web
         m[M[[2]],M[[3]]]<-web
         return(m)
      }
      return(replicate(N, sauve.fastr(web,m,M), simplify = FALSE))
   }
   else{
      sauve.fastc <- function(web,m,M) {
         web[,sample(ncol(web))]<-web
         m[M[[2]],M[[3]]]<-web
         return(m)
      }
      return(replicate(N, sauve.fastc(web,m,M), simplify = FALSE))
   }
}
