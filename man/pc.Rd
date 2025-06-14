% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pc.R
\name{pc}
\alias{pc}
\title{Participation coefficient}
\usage{
pc(network.or.subnet_mat1, subnet_mat2 = NULL, weighted = FALSE)
}
\arguments{
\item{network.or.subnet_mat1}{An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.}

\item{subnet_mat2}{A matrix representing one subnetwork.}

\item{weighted}{Logical. Default to FALSE. If TRUE, weighted degree defined as the sum of interaction strengths attached to a connector node is used.}
}
\value{
Return a numeric value representing participation coefficient for connector species.
}
\description{
Calculating participation coefficient for connector nodes.
}
\details{
In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.

The \strong{pc} function calculates the participation coefficient following Domínguez-García and Kéfi (2024). For each connector node \emph{i}, \eqn{pc_{i}} is calculated as two times the ratio between the lowest degree in both interaction subnetworks divided by the total degree of the node (\eqn{2\frac{d_{lowest}}{d_{total}}}).
Hence, the participation coefficient for all connector nodes (\eqn{PC_{c}}) is represented by the average value of all \eqn{pc_{i}}.

Two types of inputs \code{network.or.subnet_mat1} can be processed:
\itemize{
\item An "igraph" object with node attribute 'level' (0 for a-nodes, 1 for b-nodes,2 for c-nodes). If the input is a weighted network, the edge should have a 'weight' attribute.
\item Or a matrix representing subnetwork P, and must be input with \code{subnet_mat2} representing subnetwork Q.
}

If the inputs are two matrices, please make sure the rows of
\code{network.or.subnet_mat1} and \code{subnet_mat2} both represent the groups of connector species,i.e, the b-group species. If both matrices have row names, then the function matches row
names to produce connector nodes. Otherwise, row numbers are assigned to row names and matched. Within the two matrices (P and Q), columns represents a-group nodes and c-group nodes respectively.
Elements in matrices are non-zero values if two nodes are linked with or without weights, and 0 otherwise.
}
\examples{

## generate a random binary tripartite network
set.seed(12)
Net <- build_toy_net(11,15,16,0.2)
pc(Net)

##empirical network
data(PPH_Coltparkmeadow)
Net <- PPH_Coltparkmeadow
pc(Net)
set.seed(13)
library(igraph)
E(Net)$weight<-runif(length(E(Net)),0.1,1)#random weights assigned
pc(Net,weighted=T)

##input as binary matrices,with row names.
md1 <- matrix(sample(c(0,1),8*11,replace=TRUE),8,11,dimnames = list(paste0("b",1:8),paste0("c",1:11)))
md2 <- matrix(sample(c(0,1),10*12,replace=TRUE),10,12,dimnames = list(paste0("b",1:10),paste0("a",1:12)))
pc(md1,md2)

##input as weighted matrices,with row numbers as row names.
mdw1 <- matrix(sample(c(rep(0,40),runif(48,0,1))),8,11)
mdw2 <- matrix(sample(c(rep(0,40),runif(80,0,1))),10,12)
pc(mdw1,mdw2)
pc(mdw1,mdw2,weighted=T)

}
\references{
Battiston, F., Nicosia, V. & Latora, V. (2014) Structural measures for multiplex networks. Physical Review E, 89, 032804.

Domínguez-García, V. and Kéfi, S. (2024). The structure and robustness of ecological networks with two interaction types. PLOS Computational Biology, 20(1), e1011770.

Guimera, R. & Amaral, L.A.N. (2005) Cartography of complex networks: modules and universal roles. Journal of Statistical Mechanics: Theory and Experiment, 2005, P02001.
}
