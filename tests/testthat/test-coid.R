test_that("Determine the type of parameters", {


   expect_error(coid(c(1:10)),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(coid(c("a", "b", "c")),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(coid(c(T, F, F, T, NA)),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(coid(matrix(1:10,2,5),c("a", "b", "c")),
                "please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
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
   expect_error(coid(m1,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_message(coid(m1,m6),
                  "re-check whether the row name of network.or.subnet_mat1 is corresponding to the row name of subnet_mat2!!!")
   expect_error(coid(m4,m5),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(coid(m3,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
})


test_that("Input a big network data", {
   ma<-Multi_motif("all")
   for(i in 23:31){
      expect_true(class(coid(ma[[i]]))=="numeric")
   }
   set.seed(1)
   MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
   expect_error(coid(t(MAT[[3]]),t(MAT[[4]])), "missing value where TRUE/FALSE needed")

   MA<-build_net(5,3,3,0.9)
   expect_identical(class(coid(MA)), "numeric")
   set.seed(1)
   m8<-matrix(1:6,3,2)
   colnames(m8)<-paste0("species",seq=1:2)
   m9<-matrix(1:8,2,4)
   rownames(m9)<-paste0("species",seq=c(2,1))
   expect_true(is.na(coid(m8,m9)), TRUE)
})
