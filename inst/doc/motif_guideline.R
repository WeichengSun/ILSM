## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=F, out.width = "100%"-----------------------------------------------
knitr::include_graphics("../man/figure/pattern.png")

## ----echo = T,eval = F--------------------------------------------------------
#  library(ILSM)
#  type<-c("M111","M112","M113","M114","M211","M212","M213","M311","M312","M411","M121","M122-1",
#         "M122-2","M122-3","M123-1","M123-2","M123-3","M123-4","M123-5","M221-1","M221-2",
#         "M221-3","M222-1","M222-2","M222-3","M222-4","M222-5","M222-6","M222-7","M222-8",
#         "M222-9","M321-1","M321-2","M321-3","M321-4","M321-5","M131","M132-1","M132-2",
#         "M132-3","M132-4","M132-5","M231-1","M231-2","M231-3","M231-4","M231-5","M141")
#  mr <- par(mfrow=c(6,8),mar=c(1,1,3,1))
#  a<-Multi_motif("all")
#   for(i in 1:48){
#       plot(a[[i]],
#            vertex.size=30, vertex.label=NA,
#            vertex.color="#D0E7ED",main=type[i])
#  }
#  par(mr)
#  

## ----echo=F, out.width = "100%"-----------------------------------------------
knitr::include_graphics("../man/figure/motif_ILSM.png")

