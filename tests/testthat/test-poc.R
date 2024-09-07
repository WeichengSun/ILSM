test_that("Determine the type of parameters", {


   expect_error(poc(c(1:10)),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(poc(c("a", "b", "c")),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(poc(c(T, F, F, T, NA)),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(poc(matrix(1:10,2,5),c("a", "b", "c")),
                "please check the type of 'subnet_mat2' or the row numeber of this 'two matrix'")
   m1<-matrix(1:10,5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m3<-matrix(1:12,6,2)
   rownames(m3)<-c(paste0("species",seq=1:5),NA)
   m4<-matrix(1:10,5,2)
   rownames(m4)<-c(paste0("species",seq=1:4),NA)
   m5<-matrix(1:15,5,3)
   rownames(m5)<-c(paste0("species",seq=1:4),NA)
   m7<-matrix(1:18,6,3)
   rownames(m7)<-c(paste0("species",seq=c(1,3,2,5,4)),NA)
   expect_error(poc(m1,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(poc(m4,m5),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(poc(m3,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
})


test_that("Input a big network data", {
   ma<-Multi_motif("all")
   for(i in 23:31){
      expect_true(poc(ma[[i]])==1)
   }

   MAT <- build_net(11,22,21,0.2,asmatrices=TRUE)
   expect_identical(class(poc(t(MAT[[3]]),t(MAT[[4]]))), "numeric" )

   MA<-build_net(5,3,3,0.9)
   expect_identical(class(poc(MA)), "numeric")
   m8<-matrix(1:6,3,2)
   colnames(m8)<-paste0("species",seq=1:2)
   m9<-matrix(1:8,2,4)
   rownames(m9)<-paste0("species",seq=c(2,1))
   expect_identical(class(poc(m8,m9)), "numeric")
})
