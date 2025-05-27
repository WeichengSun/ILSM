Intra_guild_role <-
   function(subnet_mat1, subnet_mat2,
            AA= NULL, BB= NULL, CC= NULL,
            weighted =FALSE) {
      if (inherits(subnet_mat1, c("matrix", "data.frame")) == T &&
          inherits(subnet_mat2, c("matrix", "data.frame")) == T ){
         if(is.null(AA)){
            AA <- matrix(0, nrow = ncol(subnet_mat1), ncol = ncol(subnet_mat1))
            dimnames(AA) <- list(colnames(subnet_mat1),colnames(subnet_mat1))
         }
         if(is.null(BB)){
            BB <- matrix(0, nrow = nrow(subnet_mat1), ncol = nrow(subnet_mat1))
            dimnames(BB) <- list(rownames(subnet_mat1),rownames(subnet_mat1))
         }
         if(is.null(CC)){
            CC <- matrix(0, nrow = ncol(subnet_mat2), ncol = ncol(subnet_mat2))
            dimnames(CC) <- list(colnames(subnet_mat2),colnames(subnet_mat2))
         }

         if (inherits(AA, c("matrix", "data.frame")) == T &&
             inherits(BB, c("matrix", "data.frame")) == T &&
             inherits(CC, c("matrix", "data.frame")) == T) {
            if (is.null(rownames(subnet_mat1)) |
                is.null(rownames(subnet_mat2)) | is.null(rownames(BB)))
               stop("Make sure these matrices either have no row names or have full row names. [No NA!!!]")
            matrow <- unique(c(rownames(subnet_mat1), rownames(subnet_mat2)))

            matrow_l <- length(matrow)
            Brow <- rownames(BB)[rownames(BB) %in% matrow]
            BB_mat <- matrix(0, matrow_l, matrow_l)
            dimnames(BB_mat) <- list(matrow, matrow)
            BB_mat[Brow, Brow] <- BB[Brow, Brow]


            if (is.null(colnames(subnet_mat1)) | is.null(rownames(AA)))
               stop("Make sure these matrices either have no row and column names or have full row and column names. [No NA!!!]")
            A_l <- ncol(subnet_mat1)
            Arow <- rownames(AA)[rownames(AA) %in% colnames(subnet_mat1)]
            AA_mat <- matrix(0, A_l, A_l)
            dimnames(AA_mat) <-
               list(colnames(subnet_mat1), colnames(subnet_mat1))
            AA_mat[Arow, Arow] <- AA[Arow, Arow]

            if (is.null(colnames(subnet_mat2)) | is.null(rownames(CC)))
               stop("Make sure these matrices either have no row and column names or have full row and columnnames. [No NA!!!]")
            C_l <- ncol(subnet_mat2)
            Crow <- rownames(CC)[rownames(CC) %in% colnames(subnet_mat2)]
            CC_mat <- matrix(0, C_l, C_l)
            dimnames(CC_mat) <-
               list(colnames(subnet_mat2), colnames(subnet_mat2))
            CC_mat[Crow, Crow] <- CC[Crow, Crow]
         }
      }

      mat_1 <- matrix(0, matrow_l, A_l)
      rownames(mat_1) <- matrow
      mat_1[rownames(subnet_mat1), ] <- subnet_mat1
      mat_2 <- matrix(0, matrow_l, C_l)
      rownames(mat_2) <- matrow
      mat_2[rownames(subnet_mat2), ] <- subnet_mat2
      dimnames(mat_1) <- NULL
      dimnames(mat_2) <- NULL
      P <- t(mat_1)
      Q <- mat_2
      AA <- AA_mat
      BB <- BB_mat
      CC <- CC_mat

      PW <- P; QW <- Q; AAW <- AA; BBW <- BB; CCW <- CC
      AA <- (AA!=0)*1; BB <- (BB!=0)*1; CC <- (CC!=0)*1
      P <- (P!=0)*1; Q <- (Q!=0)*1

      AA_no <- AA == 0
      diag(AA) <- 0
      diag(AA_no) <- 0
      AA_tow <- (AA!=0)*1

      BB_no <- BB == 0
      diag(BB) <- 0
      diag(BB_no) <- 0
      BB_two <- (BB!=0)*1


      CC_no <- CC == 0
      diag(CC) <- 0
      diag(CC_no) <- 0
      CC_tow <- (CC!=0)*1


      Ob <- colSums(P); Rb <- rowSums(Q)
      Ubb <- t(P) %*% P
      Vbb <- Q %*% t(Q)
      UVbb <- Ubb * Vbb
      B2C_NO <- rowSums((Q%*%CC_no) *Q) /2
      B2C_YES <- rowSums((Q%*%CC_tow) *Q) /2
      B2A_NO <- rowSums((t(P)%*%AA_no) *t(P)) /2
      B2A_YES <- rowSums((t(P)%*%AA_tow) *t(P)) /2

      BB_role <- NULL

      M111 <- Ob * Rb
      BB_role <- cbind(BB_role,M111)

      M112_1 <- Ob * B2C_NO
      BB_role <- cbind(BB_role,M112_1)

      M112_2 <- Ob * B2C_YES
      BB_role <- cbind(BB_role,M112_2)

      M211_1 <- B2A_NO * Rb
      BB_role <- cbind(BB_role,M211_1)

      M211_2 <- B2A_YES * Rb
      BB_role <- cbind(BB_role,M211_2)

      M212_1 <- B2A_NO * B2C_NO
      BB_role <- cbind(BB_role,M212_1)

      M212_2 <- B2A_YES * B2C_NO
      BB_role <- cbind(BB_role,M212_2)

      M212_3 <- B2A_NO * B2C_YE
      BB_role <- cbind(BB_role,M212_3)

      M212_4 <- B2A_YES * B2C_YES
      BB_role <- cbind(BB_role,M212_4)

      # mat_tri <- function(mat) {sum(mat[upper.tri(mat)])}

      M121_1 <- rowSums(UVbb * BB_no)
      BB_role <- cbind(BB_role,M121_1)

      M121_2 <- rowSums(UVbb * BB_two)
      BB_role <- cbind(BB_role,M121_2)


      P_all_yes <- t(apply(P, 2, function(x) { A <- x %*% t(x); A * AA}))
      P_one_yes <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); A *AA}))
      P_two_yes <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); t(A) *AA}))

      P_all_no <- t(apply(P, 2, function(x) { A <- x %*% t(x); A * AA_no}))
      P_one_no <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); A * AA_no}))
      P_two_no <- t(apply(P, 2, function(x) { A <- (x == 0) %*% t(x); t(A) * AA_no}))

      Q_all_yes <- t(apply(Q, 1, function(x) { A <- x %*% t(x); A * CC}))
      Q_one_yes <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); A * CC}))
      Q_two_yes <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); t(A) * CC}))

      Q_all_no <- t(apply(Q, 1, function(x) { A <- x %*% t(x); A * CC_no}))
      Q_one_no <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); A * CC_no}))
      Q_two_no <- t(apply(Q, 1, function(x) { A <- (x == 0) %*% t(x); t(A) * CC_no}))

      AA_no_P_2 <- P_two_no %*% t(P_one_no)
      AA_yes_P_2 <- P_two_yes %*% t(P_one_yes)
      AA_no_P_3 <- P_all_no %*% t(P_one_no)
      AA_yes_P_3 <- P_all_yes %*% t(P_one_yes)
      AA_no_P_3t <- t(AA_no_P_3)
      AA_yes_P_3t <- t(AA_yes_P_3)
      AA_no_P_4 <- (P_all_no %*% t(P_all_no)) / 2
      AA_yes_P_4 <- (P_all_yes %*% t(P_all_yes)) / 2


      CC_no_Q_2 <- Q_two_no %*% t(Q_one_no)
      CC_yes_Q_2 <- Q_two_yes %*% t(Q_one_yes)
      CC_no_Q_3 <- Q_all_no %*% t(Q_one_no)
      CC_yes_Q_3 <- Q_all_yes %*% t(Q_one_yes)
      CC_no_Q_4 <- (Q_all_no %*% t(Q_all_no)) / 2
      CC_yes_Q_4 <- (Q_all_yes %*% t(Q_all_yes)) / 2



      M122_1_1 <- rowSums(Ubb * BB_no * CC_no_Q_2)
      BB_role <- cbind(BB_role,M122_1_1)

      M122_1_2 <- rowSums(Ubb * BB * CC_no_Q_2)
      BB_role <- cbind(BB_role,M122_1_2)

      M122_1_3 <- rowSums(Ubb * BB_no * CC_yes_Q_2)
      BB_role <- cbind(BB_role,M122_1_3)

      M122_1_4 <- rowSums(Ubb * BB * CC_yes_Q_2)
      BB_role <- cbind(BB_role,M122_1_4)

      M122_2_1 <- Ubb * BB_no * CC_no_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_1), colSums(M122_2_1))

      M122_2_2 <- Ubb * BB * CC_no_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_2), colSums(M122_2_2))

      M122_2_3 <- Ubb * BB_no * CC_yes_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_3), colSums(M122_2_3))

      M122_2_4 <- Ubb * BB * CC_yes_Q_3
      BB_role <- cbind(BB_role,rowSums(M122_2_4), colSums(M122_2_4))

      M122_3_1 <- rowSums(Ubb * BB_no * CC_no_Q_4)
      BB_role <- cbind(BB_role,M122_3_1)

      M122_3_2 <- rowSums(Ubb * BB * CC_no_Q_4)
      BB_role <- cbind(BB_role,M122_3_2)

      M122_3_3 <- rowSums(Ubb * BB_no * CC_yes_Q_4)
      BB_role <- cbind(BB_role,M122_3_3)

      M122_3_4 <- rowSums(Ubb * BB * CC_yes_Q_4)
      BB_role <- cbind(BB_role,M122_3_4)

      M221_1_1 <- rowSums(AA_no_P_2 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_1_1)

      M221_1_2 <- rowSums(AA_yes_P_2 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_1_2)

      M221_1_3 <- rowSums(AA_no_P_2 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_1_3)

      M221_1_4 <- rowSums(AA_yes_P_2 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_1_4)


      M221_2_1 <- sum(AA_no_P_3 * BB_no * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_1), colSums(M221_2_1))

      M221_2_2 <- sum(AA_yes_P_3 * BB_no * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_2), colSums(M221_2_2))

      M221_2_3 <- sum(AA_no_P_3 * BB * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_3), colSums(M221_2_3))

      M221_2_4 <- sum(AA_yes_P_3 * BB * Vbb)
      BB_role <- cbind(BB_role,rowSums(M221_2_4), colSums(M221_2_4))

      M221_3_1 <- rowSums(AA_no_P_4 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_3_1)

      M221_3_2 <- rowSums(AA_yes_P_4 * BB_no * Vbb)
      BB_role <- cbind(BB_role,M221_3_2)

      M221_3_3 <- rowSums(AA_no_P_4 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_3_3)

      M221_3_4 <- rowSums(AA_yes_P_4 * BB * Vbb)
      BB_role <- cbind(BB_role,M221_3_4)

      rownames(BB_role) <- matrow
      colnames(BB_role) <- paste0("role",1:ncol(BB_role))

      intra_guild_role <- BB_role



      if(weighted){
         All_add <- function(a,b){ return((a!=0)*(b!=0)*(a+b)) }
         One_add <- function(a,b){ return((a!=0)*(b!=0)*(a)) }

         ObW <- colSums(PW)
         RbW <- rowSums(QW)

         B2C_NOW <- apply(QW,1,function(x){ A<- apply(CC_no,2,function(y){ One_add(x,y) });  return(sum(All_add(x,t(A)))/2)})
         B2C_YESW <- apply(QW,1,function(x){ A<- apply(CCW,2,function(y){ All_add(x,y) });  return(sum(All_add(x,t(A)))/2)})

         B2A_NOW <- apply(PW,2,function(x){ A<- apply(AA_no,2,function(y){ One_add(x,y) });  return(sum(All_add(x,t(A)))/2)})
         B2A_YESW <- apply(PW,2,function(x){ A<- apply(AAW,2,function(y){ All_add(x,y) });  return(sum(All_add(x,t(A)))/2)})

         BB_role <- NULL

         M111 <- (ObW * Rb + RbW * Ob)/2
         BB_role <- cbind(BB_role,M111)

         M112_1 <- (ObW * B2C_NO + B2C_NOW * Ob  )/3
         BB_role <- cbind(BB_role,M112_1)

         M112_2 <- (ObW * B2C_YES + B2C_YESW * Ob  )/4
         BB_role <- cbind(BB_role,M112_2)

         M211_1 <- (B2A_NOW * Rb + RbW * B2A_NO )/3
         BB_role <- cbind(BB_role,M211_1)

         M211_2 <- (B2A_YESW * Rb +  RbW * B2A_YES )/4
         BB_role <- cbind(BB_role,M211_2)

         M212_1 <- (B2A_NOW * B2C_NO + B2C_NOW *B2A_NO)/4
         BB_role <- cbind(BB_role,M212_1)

         M212_2 <- (B2A_YESW * B2C_NO + B2C_NOW * B2A_YES)/5
         BB_role <- cbind(BB_role,M212_2)

         M212_3 <- (B2A_NOW * B2C_YES + B2C_YESW * B2A_NO)/5
         BB_role <- cbind(BB_role,M212_3)

         M212_4 <- (B2A_YESW * B2C_YES + B2C_YESW * B2A_YES)/6
         BB_role <- cbind(BB_role,M212_4)


         A_2BW <- t(apply(PW,2,function(x){apply(PW,2,function(y){ sum(All_add(x,y)) })}))
         C_2BW <- t(apply(QW,1,function(x){apply(QW,1,function(y){ sum(All_add(x,y)) })}))


         M121_1 <- rowSums( (A_2BW*Vbb + C_2BW*Ubb) * BB_no)/4
         BB_role <- cbind(BB_role,M121_1)

         M121_2 <- rowSums((A_2BW*Vbb + C_2BW*Ubb)* BB_two + UVbb * BBW )/5
         BB_role <- cbind(BB_role,M121_2)


         PW_all_yes <- t(apply(PW, 2, function(x){A <- (x!=0)%*%t(x); All_add(A,x)}))  # All_add(t(All_add(t(All_add(A,x)),x)),PP)
         PW_one_yes <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); A }))
         PW_two_yes <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); t(A) }))

         PW_all_no <- t(apply(PW, 2, function(x){A <- (x!=0)%*%t(x);  All_add(A,x) * AA_no}))  # All_add(t(All_add(t(All_add(A,x)),x)),PP)
         PW_one_no <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); A * AA_no}))
         PW_two_no <- t(apply(PW, 2, function(x){A <- (x==0)%*%t(x); t(A) * AA_no}))

         QW_all_yes <- t(apply(QW, 1, function(x){A <- (x!=0)%*%t(x); All_add(A,x)}))
         QW_one_yes <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); A}))
         QW_two_yes <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); t(A) }))

         QW_all_no <- t(apply(QW, 1, function(x){A <- (x!=0)%*%t(x); All_add(A,x)* CC_no}))
         QW_one_no <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); A * CC_no}))
         QW_two_no <- t(apply(QW, 1, function(x){A <- (x==0)%*%t(x); t(A) * CC_no}))




         AA_no_PW_2 <- t(apply(PW_two_no,1,function(x){apply(PW_one_no,1,function(y){ sum(All_add(x,y))})}))
         AA_yes_PW_2 <- t(apply(PW_two_yes,1,function(x){apply(PW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_tow)*(x+y+as.vector(AAW)))})}))

         AA_no_PW_3 <- t(apply(PW_all_no,1,function(x){apply(PW_one_no,1,function(y){ sum(All_add(x,y))})}))
         AA_yes_PW_3 <- t(apply(PW_all_yes,1,function(x){apply(PW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_tow)*(x+y+as.vector(AAW)))})}))

         AA_no_PW_3t <- t(AA_no_PW_3)
         AA_yes_PW_3t <- t(AA_yes_PW_3)

         AA_no_PW_4 <- t(apply(PW_all_no,1,function(x){apply(PW_all_no,1,function(y){ sum(All_add(x,y))/2})}))
         AA_yes_PW_4 <- t(apply(PW_all_yes,1,function(x){apply(PW_all_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(AA_tow)*(x+y+as.vector(AAW))) /2})}))


         CC_no_QW_2 <- t(apply(QW_two_no,1,function(x){apply(QW_one_no,1,function(y){ sum(All_add(x,y))})}))
         CC_yes_QW_2 <- t(apply(QW_two_yes,1,function(x){apply(QW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_tow)*(x+y+as.vector(CCW)))})}))

         CC_no_QW_3 <- t(apply(QW_all_no,1,function(x){apply(QW_one_no,1,function(y){ sum(All_add(x,y))})}))
         CC_yes_QW_3 <- t(apply(QW_all_yes,1,function(x){apply(QW_one_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_tow)*(x+y+as.vector(CCW)))})}))

         CC_no_QW_4 <- t(apply(QW_all_no,1,function(x){apply(QW_all_no,1,function(y){ sum(All_add(x,y))/2 })}))
         CC_yes_QW_4 <- t(apply(QW_all_yes,1,function(x){apply(QW_all_yes,1,function(y){ sum((x!=0)*(y!=0)*as.vector(CC_tow)*(x+y+as.vector(CCW)))/2 })}))





         M122_1_1 <- rowSums((A_2BW * CC_no_Q_2 + CC_no_QW_2 * Ubb) * BB_no)/4
         BB_role <- cbind(BB_role,M122_1_1)

         M122_1_2 <- rowSums((A_2BW * CC_no_Q_2 + CC_no_QW_2 * Ubb) * BB_two + Ubb * CC_no_Q_2 * BBW)/5
         BB_role <- cbind(BB_role,M122_1_2)

         M122_1_3 <- rowSums((A_2BW * CC_yes_Q_2 + CC_yes_QW_2 * Ubb) * BB_no)/5
         BB_role <- cbind(BB_role,M122_1_3)

         M122_1_4 <- rowSums((A_2BW * CC_yes_Q_2 + CC_yes_QW_2 * Ubb) * BB_two + Ubb * CC_yes_Q_2 * BBW)/6
         BB_role <- cbind(BB_role,M122_1_4)

         M122_2_1 <- ((A_2BW * CC_no_Q_3 + CC_no_QW_3 * Ubb) * BB_no)/5
         BB_role <- cbind(BB_role,rowSums(M122_2_1), colSums(M122_2_1))

         M122_2_2 <- ((A_2BW * CC_no_Q_3 + CC_no_QW_3 * Ubb) * BB_two + Ubb * CC_no_Q_3 * BBW)/6
         BB_role <- cbind(BB_role,rowSums(M122_2_2), colSums(M122_2_2))

         M122_2_3 <- ((A_2BW * CC_yes_Q_3 + CC_yes_QW_3 * Ubb) * BB_no)/6
         BB_role <- cbind(BB_role,rowSums(M122_2_3), colSums(M122_2_3))

         M122_2_4 <- ((A_2BW * CC_yes_Q_3 + CC_yes_QW_3 * Ubb) * BB_two + Ubb * CC_yes_Q_3 * BBW)/7
         BB_role <- cbind(BB_role,rowSums(M122_2_4), colSums(M122_2_4))



         M122_3_1 <- rowSums((A_2BW * CC_no_Q_4 + CC_no_QW_4 * Ubb) * BB_no)/6
         BB_role <- cbind(BB_role,M122_3_1)

         M122_3_2 <- rowSums((A_2BW * CC_no_Q_4 + CC_no_QW_4 * Ubb) * BB_two + Ubb * CC_no_Q_4 * BBW)/7
         BB_role <- cbind(BB_role,M122_3_2)

         M122_3_3 <- rowSums((A_2BW * CC_yes_Q_4 + CC_yes_QW_4 * Ubb) * BB_no)/7
         BB_role <- cbind(BB_role,M122_3_3)

         M122_3_4 <- rowSums((A_2BW * CC_yes_Q_4 + CC_yes_QW_4 * Ubb) * BB_two + Ubb * CC_yes_Q_4 * BBW)/8
         BB_role <- cbind(BB_role,M122_3_4)


         M221_1_1 <- rowSums((AA_no_PW_2 * Vbb + C_2BW * AA_no_P_2) * BB_no)/4
         BB_role <- cbind(BB_role,M221_1_1)

         M221_1_2 <- rowSums((AA_yes_PW_2 * Vbb + C_2BW * AA_yes_P_2) * BB_no)/5
         BB_role <- cbind(BB_role,M221_1_2)

         M221_1_3 <- rowSums((AA_no_PW_2 * Vbb + C_2BW * AA_no_P_2) * BB_two + AA_no_P_2 * Vbb * BBW)/5
         BB_role <- cbind(BB_role,M221_1_3)

         M221_1_4 <- rowSums((AA_yes_PW_2 * Vbb + C_2BW * AA_yes_P_2) * BB_two + AA_yes_P_2 * Vbb * BBW)/6
         BB_role <- cbind(BB_role,M221_1_4)



         M221_2_1 <- ((AA_no_PW_3 * Vbb + C_2BW * AA_no_P_3) * BB_no)/5
         BB_role <- cbind(BB_role,rowSums(M221_2_1), colSums(M221_2_1))

         M221_2_2 <- ((AA_yes_PW_3 * Vbb + C_2BW * AA_yes_P_3) * BB_no)/6
         BB_role <- cbind(BB_role,rowSums(M221_2_2), colSums(M221_2_2))

         M221_2_3 <- ((AA_no_PW_3 * Vbb + C_2BW * AA_no_P_3) * BB_two + AA_no_P_3 * Vbb * BBW)/6
         BB_role <- cbind(BB_role,rowSums(M221_2_3), colSums(M221_2_3))

         M221_2_4 <- ((AA_yes_PW_3 * Vbb + C_2BW * AA_yes_P_3) * BB_two + AA_yes_P_3 * Vbb * BBW)/7
         BB_role <- cbind(BB_role,rowSums(M221_2_4), colSums(M221_2_4))


         M221_3_1 <- rowSums((AA_no_PW_4 * Vbb + C_2BW * AA_no_P_4) * BB_no)/6
         BB_role <- cbind(BB_role,M221_3_1)

         M221_3_2 <- rowSums((AA_yes_PW_4 * Vbb + C_2BW * AA_yes_P_4) * BB_no)/7
         BB_role <- cbind(BB_role,M221_3_2)

         M221_3_3 <- rowSums((AA_no_PW_4 * Vbb + C_2BW * AA_no_P_4) * BB_two + AA_no_P_4 * Vbb * BBW)/7
         BB_role <- cbind(BB_role,M221_3_3)

         M221_3_4 <- rowSums((AA_yes_PW_4 * Vbb + C_2BW * AA_yes_P_4) * BB_two + AA_yes_P_4 * Vbb * BBW)/8
         BB_role <- cbind(BB_role,M221_3_4)

         rownames(BB_role) <- matrow
         colnames(BB_role) <- paste0("weighted role",1:ncol(BB_role))
         intra_guild_mean_role <- BB_role

         intra_guild_mean_role <- intra_guild_mean_role/intra_guild_role

         intra_guild_mean_role <- replace(intra_guild_mean_role, which(is.nan(intra_guild_mean_role)),NA)
         return(intra_guild_mean_role)
      }
      else
         return(intra_guild_role)
   }
