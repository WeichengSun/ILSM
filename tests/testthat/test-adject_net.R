test_that("Determine whether a network can be adapted into an interconnected network", {
   m1<-matrix(c(0,0,0,1,0,0),3,2)
   m2<-matrix(c(0,0,0,2,0,0),2,3)
   N<-igraph_from_matrices(m1,m2)
   expect_error(icmotif_count(N),
                "'network' is not a 'interconnecting' network by central layer")

})
