test_that("Determine the type of parameters", {
#' @srrstats {G5.8,G5.8b,G5.8c,G5.8d} The test data include: numeric vector, 'NA' data and data with less fields(less than minimum).

   expect_error(icmotif_count(c(1:10)),
                "Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   expect_error(icmotif_count(c("a", "b", "c")),
                "Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   expect_error(icmotif_count(c(T, F, F, T, NA)),
                "Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   expect_error(icmotif_count(matrix(1:10,2,5),c("a", "b", "c")),
                "Error: please check the tyep of network.or.subnet_mat1 and other parameters!!!")
   m1<-matrix(1:10,5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m3<-matrix(1:12,6,2)
   rownames(m3)<-c(paste0("species",seq=1:5),NA)
   m4<-matrix(1:10,5,2)
   rownames(m4)<-c(paste0("species",seq=1:4),NA)
   m5<-matrix(1:15,5,3)
   rownames(m5)<-c(paste0("species",seq=1:4),NA)
   m6<-matrix(1:18,6,3)
   rownames(m6)<-paste0("species",seq=2:7)
   m7<-matrix(1:18,6,3)
   rownames(m7)<-c(paste0("species",seq=c(1,3,2,5,4)),NA)
   expect_error(icmotif_count(m1,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_message(icmotif_count(m1,m6),
                  "re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
   expect_error(icmotif_count(m4,m5),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(icmotif_count(m3,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
})


test_that("Input a big network data", {
   ma<-Multi_motif("all")
   for(i in 23:31){
      expect_error(icmotif_count(ma[[i]]),
                   "Error: please input a large 'number of interconnecting species >=4' network data!!!")
   }

   MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
   expect_error(icmotif_count(t(MAT[[3]]),t(MAT[[4]])),
                "Error: please input a large 'number of interconnecting species >=4' network data!!!")

   MA<-build_net(5,3,3,0.9)
   expect_error(icmotif_count(MA),
                "Error: please input a large 'number of interconnecting species >=4' network data!!!")
   m8<-matrix(1:6,3,2)
   colnames(m8)<-paste0("species",seq=1:2)
   m9<-matrix(1:8,2,4)
   rownames(m9)<-paste0("species",seq=c(2,1))
   expect_error(icmotif_count(m8,m9),
                "Error: please input a large 'number of interconnecting species >=4' network data!!!")
})
