test_that("Determine the type of matrices of inputing", {
   m1<-matrix(1:10,5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(1:12,4,3)
   m3<-matrix(1:12,6,2)
   rownames(m3)<-c(paste0("species",seq=1:5),NA)
   m4<-matrix(1:10,5,2)
   rownames(m4)<-c(paste0("species",seq=1:4),NA)
   m5<-matrix(1:15,5,3)
   rownames(m5)<-c(paste0("species",seq=1:4),NA)
   m6<-matrix(1:15,5,3)
   rownames(m6)<-paste0("species",seq=1:5)
   m7<-matrix(1:18,6,3)
   rownames(m7)<-c(paste0("species",seq=c(1,3,2,5,4)),NA)
   m8<-matrix(1:18,6,3)
   rownames(m8)<-c(paste0("species",seq=c(1,3,2,4)),NA,NA)
   expect_message(igraph_from_matrices(m1,m2),
                  "You will obtain a network that contains non-connecting node!!!")
   expect_message(igraph_from_matrices(m1,m6),
                  "Because you enter two matrices with the same number of rows, the function will default to them having the same row names. Informed!!")
   expect_error(igraph_from_matrices(m1,m7),
                "Make sure matrices either have no row names or have full row names. No NA!!!")
   expect_error(igraph_from_matrices(m4,m5),
                "Error: There is NA in the column name of mat1 or the row name of mat2!!!")
   expect_error(igraph_from_matrices(m3,m7),
                "Error: There is NA in the column name of mat1 or the row name of mat2!!!")
   expect_error(igraph_from_matrices(m3,m8),
                "Error: please check whether the column name of mat1 is corresponding to the row name of mat2!!!")
})


test_that("Make sure the function is implemented", {
   m1<-matrix(sample(c(rep(1,7),rep(0,3))),5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(sample(c(rep(1,9),rep(0,6))),5,3)
   rownames(m2)<-c(paste0("species",seq=c(1,3,2,5,4)))
   m3<-matrix(sample(c(rep(1,8),rep(0,4))),6,2)
   m4<-matrix(sample(c(rep(1,9),rep(0,6))),5,3)
   N<-igraph_from_matrices(m1,m2)
   M<-igraph_from_matrices(m3,m4)
   m<-igraph_from_matrices(m4,m3)
   expect_identical(class(N),
                    "igraph")
   expect_identical(length(N),
                    10L)
   expect_identical(length(M),
                    11L)
   expect_true(length(M)==length(m))
   expect_length(motif_count(N),
                    44L)
   expect_true(motif_count(M)[1]==motif_count(m)[1])
})
