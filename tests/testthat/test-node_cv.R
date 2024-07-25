test_that("Determine the type of matrices of inputing", {

   expect_error(node_cv(c(1:10)),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(node_cv(c("a", "b", "c")),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(node_cv(c(T, F, F, T, NA)),
                "please check the type of 'network.or.subnet_mat1'")
   expect_error(node_cv(matrix(1:10,2,5),c("a", "b", "c")),
                "please check the type of 'network.or.subnet_mat1'")

   m1<-matrix(1:10,5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m3<-matrix(1:12,6,2)
   rownames(m3)<-c(paste0("species",seq=1:5),NA)
   m4<-matrix(1:10,5,2)
   rownames(m4)<-c(paste0("species",seq=1:4),NA)
   m5<-matrix(1:15,5,3)
   rownames(m5)<-c(paste0("species",seq=1:4),NA)
   # m6<-matrix(1:18,6,3)
   # rownames(m6)<-paste0("species",seq=3:8)
   m7<-matrix(1:18,6,3)
   rownames(m7)<-c(paste0("species",seq=c(1,3,2,5,4)),NA)
   m8<-matrix(1:18,6,3)
   rownames(m8)<-c(paste0("species",seq=c(1,3,2,4)),NA,NA)
   expect_error(node_cv(m1,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(node_cv(m4,m5),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(node_cv(m3,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(node_cv(m3,m8),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(node_cv(m8,m3),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
})


test_that("Make sure the function is implemented", {
   m1<-matrix(sample(c(rep(1,5),rep(0,5))),5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(sample(c(rep(1,7),rep(0,8))),5,3)
   rownames(m2)<-c(paste0("species",seq=c(1,3,2,5,4)))
   N<-node_cv(m1,m2)
   expect_identical(class(N),
                    "data.frame")
   expect_identical(ncol(N),
                    8L)
   expect_length(rownames(N),
                 10L)
})
