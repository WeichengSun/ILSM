
test_that("Determine the type of result", {
   type=c("M111","M112","M113","M211","M212","M213","M311","M312","M121-1","M122-1","M122-2",
          "M122-3","M123-1","M123-2","M123-3","M123-4","M123-5","M221-1","M221-2","M221-3","M222-1","M222-2",
          "M222-3","M222-4","M222-5","M222-6","M222-7","M222-8","M321-1","M321-2","M321-3","M321-4","M321-5",
          "M131","M132-1","M132-2","M132-3","M132-4","M132-5","M231-1","M231-2","M231-3","M231-4","M231-5")
   for(i in 1:44){
      ms<-Multi_motif(type[i])
      expect_true(inherits(ms,"igraph"))
   }
   ma<-Multi_motif("all")
   expect_true(inherits(ma,"list"))
   for(i in 1:44)
      expect_true(inherits(ma[[i]],"igraph"))
})


test_that("Determine the correction of single motif and all motifs", {
   type=c("M111","M112","M113","M211","M212","M213","M311","M312","M121-1","M122-1","M122-2",
          "M122-3","M123-1","M123-2","M123-3","M123-4","M123-5","M221-1","M221-2","M221-3","M222-1","M222-2",
          "M222-3","M222-4","M222-5","M222-6","M222-7","M222-8","M321-1","M321-2","M321-3","M321-4","M321-5",
          "M131","M132-1","M132-2","M132-3","M132-4","M132-5","M231-1","M231-2","M231-3","M231-4","M231-5")
   ma<-Multi_motif("all")
   for(i in 1:44){
      ms<-Multi_motif(type[i])
      expect_false(sum(!(as.matrix(ma[[i]])==as.matrix(ms)))!=0)
   }
})
