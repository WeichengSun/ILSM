test_that("Build a large network", {
  expect_error(build_net(4,2,3,0.5),
               "Error: please make lay_0>=3, lay_1>=3 and lay_2>=3!!!")
   expect_error(build_net(4,2,3,0.5,asmatrices=TRUE),
                "Error: please make lay_0>=3, lay_1>=3 and lay_2>=3!!!")
   expect_error(build_net(2,4,3,0.5),
                "Error: please make lay_0>=3, lay_1>=3 and lay_2>=3!!!")
   expect_error(build_net(4,3,2,0.5),
                "Error: please make lay_0>=3, lay_1>=3 and lay_2>=3!!!")
   expect_error(build_net(2,2,2,0.5),
                "Error: please make lay_0>=3, lay_1>=3 and lay_2>=3!!!")
   expect_equal(length(build_net(5,5,5,0.5)),
                15)
   N<-build_net(5,5,5,0.5,TRUE)
   expect_named(N,c( "network","supraadjacency_matrix","subnetwork1","subnetwork2" ))
})

test_that("Ensure the accuracy of the conversion into a matrix", {
   set.seed(1)
   M<-build_net(5,5,5,0.5)
   set.seed(1)
   N<-build_net(5,5,5,0.5,TRUE)
   expect_identical(class(M),"igraph")
   expect_identical(class(N[[1]]),"igraph")
   expect_identical(class(N[[2]]),c("matrix","array"))
   expect_identical(class(N[[3]]),c("matrix","array"))
   expect_identical(class(N[[4]]),c("matrix","array"))
   expect_true(sum(!(as.matrix(M)==as.matrix(N[[1]])))==0)
   expect_identical(sum(!(as.matrix(N[[1]])==N[[2]]))==0,TRUE)
   expect_identical(sum(rownames(N[[3]])==rownames(N[[4]])),5L)
   expect_identical(sum(sort(rownames(N[[2]]))==sort(c(rownames(N[[3]]),colnames(N[[3]]),colnames(N[[4]])))),15L)
})

