SavueR<-function(N,m,type=c("row","col")){
   M<-emptyR(m)
   web <- as.matrix(M[[1]])
   if(type=="row"){
      savue.fastr <- function(web,m,M) {
         web[sample(sample(nrow(web))),]<-web
         m[M[[2]],M[[3]]]<-web
         return(m)
      }
      return(replicate(N, savue.fastr(web,m,M), simplify = FALSE))
   }
   else{
      savue.fastc <- function(web,m,M) {
         web[,sample(ncol(web))]<-web
         m[M[[2]],M[[3]]]<-web
         return(m)
      }
      return(replicate(N, savue.fastc(web,m,M), simplify = FALSE))
   }
}
