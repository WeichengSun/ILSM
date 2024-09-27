test_that("Determine the type of null model", {

   #' @srrstats {G5.5} This test runs with a fixed random seed.
   ma<-Multi_motif("all")
   for(i in 23:31){
      expect_error(null_model(ma[[i]],null_type = "a"),
                   "Error: null_type is not a valid input!")
   }
   set.seed(1)
   MA<-build_net(5,3,3,0.9)
   expect_true(class(null_model(MA)[[1]])=="igraph")

})





test_that("Make sure the function is implemented", {
   m1<-matrix(sample(c(rep(1,9),rep(0,1))),5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(sample(c(rep(1,13),rep(0,2))),5,3)
   rownames(m2)<-c(paste0("species",seq=c(1,3,2,5,4)))
   N<-igraph_from_matrices(m1,m2)
   expect_identical(class(N),
                    c("igraph"))
   expect_identical(length(null_model(N,number = 4)),
                    4L)
   expect_length(null_model(N,number = 5),
                 5L)
})
