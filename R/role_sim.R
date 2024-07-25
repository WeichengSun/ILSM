#' Similarity of roles of interconnecting species
#'
#' The average of Similarity of 65 roles of interconnecting species of multilayer network.
#'
#' @param spe_role_mat A matrix of 65 columns representing the roles of interconnecting species in the motifs. Columns names are Role1, Role2, Role3 ... Role70.
#'
#' @details
#'
#' \strong{spe_role_mat}
#'
#' Should acquire from function \code{motif_role}.
#'
#' @return
#' Return a numeric value.
#'
#' @export
#'
#' @examples
#'
#' set.seed(12)
#' d <- build_net(11,22,21,0.2)
#' mr <- icmotif_role(d)
#' role_sim(mr)
#'
#' set.seed(1)
#' D <- build_net(11,22,21,0.2)
#' role_sim(icmotif_role(D))
#'
role_sim<-function(spe_role_mat){
  logi<-(rowSums(spe_role_mat)!=0)
  spe_role_mat<-spe_role_mat[logi,]
  role_sim<-NULL
  for(i in 1:(nrow(spe_role_mat)-1)){
    for(j in (i+1):nrow(spe_role_mat)){
      mat<-spe_role_mat[c(i,j),]
      role_sim<-c(role_sim,sum(apply(mat,2,min))/sum(apply(mat,2,max)))
    }
  }
  return(mean(role_sim))
}
