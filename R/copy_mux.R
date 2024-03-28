#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix drop0
#' @importFrom Matrix Diagonal
#' @importFrom Matrix Matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph distances
#'

sumR_r<-function (A, n)
{
   if (n == 1) {
      return(rbind(Matrix::colSums(A)))
   }
   else if (n == 2) {
      return(cbind(Matrix::rowSums(A)))
   }
   else {
      stop("ERROR! Not a valid dimension.")
   }
}

diagR_r<- function(x, n, offset = 0) {
   # cannot call this function just "diag", because it's an R builtin command
   # (to-do) add check on the length of x for offset != 0
   x <- as.numeric(x)
   M <- matrix(0, n, n)
   M[which(row(M) + offset == col(M))] <- x
   return(M)
}

GetLargestEigenv_r<-function (Matrix)
{
   tmp <- eigen(Matrix, only.values = F)
   lm <- which.max(Re(tmp$values[abs(Im(tmp$values)) < 1e-06]))
   tmp$vectors <- tmp$vectors[, lm]
   tmp$values <- tmp$values[lm]
   if (all(Im(tmp$vectors) == 0)) {
      tmp$vectors <- Re(tmp$vectors)
   }
   if (all(Im(tmp$values) == 0)) {
      tmp$values <- Re(tmp$values)
   }
   tmp$vectors[Matrix::which(tmp$vectors > -1e-12 & tmp$vectors <
                                1e-12)] <- 0
   if (all(tmp$vectors[Matrix::which(tmp$vectors != 0)] < 0)) {
      tmp$vectors <- -tmp$vectors
   }
   return(list(QMatrix = tmp$vectors, LMatrix = tmp$values))
}

BuildSupraAdjacencyMatrixFromExtendedEdgelist_r<-function (mEdges, Layers, Nodes, isDirected)
{
   if (max(max(mEdges[, 2]), max(mEdges[, 4])) != Layers) {
      stop("Error: expected number of layers does not match the data. Aborting process.")
   }
   edges <- data.frame(from = mEdges[, 1] + Nodes * (mEdges[,
                       2] - 1), to = mEdges[, 3] + Nodes * (mEdges[, 4] - 1),
                       weight = mEdges[, 5])
   M <- Matrix::sparseMatrix(i = edges$from, j = edges$to,
                             x = edges$weight, dims = c(Nodes * Layers, Nodes * Layers))
   if (sum(abs(M - Matrix::t(M))) > 1e-12 && isDirected == FALSE)
      M <- (M + Matrix::t(M))/2
   return(M)
}

GetMultiPageRankCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   return(GetMultiRWCentrality_r(SupraAdjacencyMatrix, Layers,
                                 Nodes, Type = "pagerank"))
}


GetMultiRWCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes, Type = "pagerank",
                                  Method = "multilayer")
{
   SupraTransitionMatrix <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix_r(SupraAdjacencyMatrix,
                                                                                 Layers, Nodes, Type = Type)
   tmp <- GetLargestEigenv_r(Matrix::t(SupraTransitionMatrix))
   LeadingEigenvector <- tmp$QMatrix
   LeadingEigenvalue <- tmp$LMatrix
   if (abs(LeadingEigenvalue - 1) > 1e-06) {
      stop(paste("  GetRWOverallOccupationProbability: ERROR! Expected leading eigenvalue equal to 1, obtained",
                 LeadingEigenvalue, ". Aborting process."))
   }
   CentralityVector <- LeadingEigenvector/sum(LeadingEigenvector)
   if (Method == "multilayer") {
      CentralityVector <- sumR_r(matrix(CentralityVector,
                                        Nodes, Layers), 2)
   }
   CentralityVector <- CentralityVector/max(CentralityVector)
   return(CentralityVector)
}


BuildSupraTransitionMatrixFromSupraAdjacencyMatrix_r<-function (SupraAdjacencyMatrix, Layers, Nodes, Type = "pagerank", r = NULL) {
   Type <- tolower(Type)
   Order <- Layers * Nodes
   SupraStrengthMatrix <- sumR_r(SupraAdjacencyMatrix, 2)
   DisconnectedNodes <- length(SupraStrengthMatrix[SupraStrengthMatrix ==
                                                      0])
   in_layer_disconnected_count <- 0
   if (Type == "pagerank" || is.null(Type) ) {
      SupraStrengthMatrix[SupraStrengthMatrix > 0] <- 1/SupraStrengthMatrix[SupraStrengthMatrix >
                                                                               0]
      SupraStrengthMatrix <- diagR_r(SupraStrengthMatrix, Order,
                                     0)
      SupraTransitionMatrix <- SupraStrengthMatrix %*% SupraAdjacencyMatrix
   }

   alpha <- NA
   if (Type == "pagerank")
      alpha <- 0.85
   else
      alpha <- 1
   if (DisconnectedNodes > 0) {
      DisconnectedNodes.ids <- which(sumR_r(SupraTransitionMatrix,
                                            2) == 0)
      if (Type == "pagerank") {
         SupraTransitionMatrix[DisconnectedNodes.ids, ] <- 1/Order
         SupraTransitionMatrix[-DisconnectedNodes.ids, ] <- alpha *
            SupraTransitionMatrix[-DisconnectedNodes.ids,
            ] + (1 - alpha)/Order
      }
   }
   DisconnectedNodes <- 0
   if (abs(sumR_r(sumR_r(SupraTransitionMatrix, 2), 1) - Order +
           DisconnectedNodes) > 1e-06) {
      stop(paste("  BuildSupraTransitionMatrixFromSupraAdjacencyMatrix: ERROR! Problems in building the supra-transition matrix -> ",
                 abs(sumR_r(sumR_r(SupraTransitionMatrix, 2), 1) - Order +
                        DisconnectedNodes), ". Aborting process."))
   }
   return(Matrix::drop0(SupraTransitionMatrix))
}

GetMultiHubCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   SupraMatrix <- SupraAdjacencyMatrix %*% Matrix::t(SupraAdjacencyMatrix)
   tmp <- GetLargestEigenv_r(SupraMatrix + 1e-16)
   LeadingEigenvector <- tmp$QMatrix
   CentralityVector <- sumR_r(matrix(LeadingEigenvector, Nodes,
                                     Layers), 2)
   CentralityVector <- CentralityVector/max(CentralityVector)
   return(CentralityVector)
}

GetMultiAuthCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   SupraMatrix <- Matrix::t(SupraAdjacencyMatrix) %*% SupraAdjacencyMatrix
   tmp <- GetLargestEigenv_r(SupraMatrix + 1e-16)
   LeadingEigenvector <- tmp$QMatrix
   CentralityVector <- sumR_r(matrix(LeadingEigenvector, Nodes,
                                     Layers), 2)
   CentralityVector <- CentralityVector/max(CentralityVector)
   return(CentralityVector)
}

GetMultiKatzCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   tmp <- GetLargestEigenv_r(Matrix::t(SupraAdjacencyMatrix))
   LeadingEigenvalue <- tmp$LMatrix
   deltaTensor <- kronecker(Diagonal(Nodes,1), Diagonal(Layers,1))
   a <- 0.99999/abs(LeadingEigenvalue)
   KatzKernelTensor <- solve(deltaTensor - a * SupraAdjacencyMatrix)
   KatzCentralitySupraVector <- KatzKernelTensor %*% Matrix(1,Nodes *
                                                               Layers, 1)
   CentralityVector <- sumR_r(matrix(KatzCentralitySupraVector,
                                     Nodes, Layers), 2)
   CentralityVector <- CentralityVector/max(CentralityVector)
   return(CentralityVector)
}


GetMultiEigenvectorCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   tmp <- GetLargestEigenv_r(Matrix::t(SupraAdjacencyMatrix))
   LeadingEigenvector <- tmp$QMatrix
   CentralityVector <- sumR_r(matrix(LeadingEigenvector, Nodes,
                                     Layers), 2)
   CentralityVector <- CentralityVector/max(CentralityVector)
   return(CentralityVector)
}

GetMultiClosenessCentrality_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   CentralityVector <- GetMultiPathStatistics_r(SupraAdjacencyMatrix,
                                                Layers, Nodes)["closeness"]
   return(CentralityVector)
}

GetMultiPathStatistics_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   if (Layers == 1) {
      g <- igraph::graph_from_adjacency_matrix(SupraAdjacencyMatrix,
                                               weighted = T)
      DM <- igraph::distances(g, mode = "all")
   }
   else {
      g.ext <- igraph::graph_from_adjacency_matrix(SupraAdjacencyMatrix,
                                                   weighted = T)
      DM <- igraph::distances(g.ext, mode = "all")
   }
   if (Layers > 1) {
      DM.blocks <- SupraAdjacencyToBlockTensor_r(DM, Layers,
                                                 Nodes)
      DM.min <- DM.blocks[[1, 1]]
      for (l1 in 1:Layers) {
         for (l2 in l1:Layers) {
            DM.min <- pmin(DM.min, DM.blocks[[l1, l2]])
         }
      }
      DM <- DM.min
   }
   closeness <- unlist(lapply(1:Nodes, function(n)
      mean(1/DM[n, ][-n])))
   avg.path.length <- 1/mean(closeness)
   return(list(distance.matrix = DM, avg.path.length = avg.path.length,
               closeness = closeness))
}


SupraAdjacencyToBlockTensor_r<-function (SupraAdjacencyMatrix, Layers, Nodes)
{
   BlockTensor <- matrix(list(), Layers, Layers)
   lapply(1:Layers, function(i) {
      lapply(1:Layers, function(j) {
         BlockTensor[[i, j]] <<- SupraAdjacencyMatrix[(1 + (i - 1) * Nodes):(i * Nodes), (1 + (j - 1) * Nodes):(j * Nodes)]
      })
   })
   return(BlockTensor)
}
