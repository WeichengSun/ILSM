

emptyR<-function(web){
   web[is.na(web)] <- 0
   if (NCOL(web) == 1 | NROW(web) == 1) {
      if (NCOL(web) == 1 & NROW(web) != 1) {
         nr <- sum(web > 0)
         nc <- 1
      }
      if (NROW(web) == 1 & NCOL(web) != 1) {
         nc <- sum(web > 0)
         nr <- 1
      }
      if (NROW(web) == 1 & NCOL(web) == 1) {
         nr <- 1
         nc <- 1
      }
   }
   cempty <- which(colSums(web) == 0)
   rempty <- which(rowSums(web) == 0)
   cind <- if (length(cempty) == 0)
      1:NCOL(web)
   else (1:NCOL(web))[-cempty]
   rind <- if (length(rempty) == 0)
      1:NROW(web)
   else (1:NROW(web))[-rempty]
   out <- web[rind, cind, drop = FALSE]
   return(list(out,rind,cind))
}


vaznullR<-function(N,m){
   M<-emptyR(m)
   web <- as.matrix(M[[1]])
   vaznull.fast <- function(web,m,M) {
      rs.p <- rowSums(web)/sum(web)
      cs.p <- colSums(web)/sum(web)
      P <- P1 <- tcrossprod(rs.p, cs.p)
      finalmat <- matrix(0, nrow(web), ncol(web))
      n.int.finalmat <- 0
      while (n.int.finalmat < sum(dim(web))) {
         sel <- sample(1:length(web), 1, prob = P, replace = TRUE)
         finalmat[sel] <- 1
         P[outer(1 * (rowSums(finalmat) > 0), 1 * (colSums(finalmat) >
                                                      0)) == 1] <- 0
         n.int.finalmat <- sum(rowSums(finalmat) > 0) + sum(colSums(finalmat) >
                                                               0)
      }
      conn.remain <- sum(web > 0) - sum(finalmat > 0)
      if (conn.remain > 0) {
         if (length(which(finalmat == 0)) == 1) {
            add <- which(finalmat == 0)
         }
         else {
            add <- sample(which(finalmat == 0), conn.remain,
                          prob = P1[finalmat == 0])
         }
         finalmat[add] <- 1
      }
      int.remain <- sum(web) - sum(finalmat)
      if (int.remain > 0) {
         add <- sample(which(finalmat > 0), int.remain, prob = P1[which(finalmat >
                                                                           0)], replace = TRUE)
         finalmat[as.numeric(names(table(add)))] <- finalmat[as.numeric(names(table(add)))] +
            table(add)
      }
      m[M[[2]],M[[3]]]<-finalmat
      return(m)
   }
   replicate(N, vaznull.fast(web,m,M), simplify = FALSE)
}
