% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hc.R
\name{hc}
\alias{hc}
\title{The proportion of hub connectors in the shared sets of nodes}
\usage{
hc(network.or.subnet_mat1, subnet_mat2 = NULL, hubs = 0.2, weighted = FALSE)
}
\arguments{
\item{network.or.subnet_mat1}{An igraph object or matrix. An "igraph" object with node attribute 'level' or a matrix representing one subnetwork. See details.}

\item{subnet_mat2}{A matrix representing one subnetwork.}

\item{hubs}{The top x\% of the nodes in the shared set of nodes with the highest degree. default to 20\%.}

\item{weighted}{Logical. Default to FALSE. If TRUE, weighted degree defined as the sum of link strengths attached to a connector node is used.}
}
\value{
Return a numeric value representing the proportion of connector nodes in node hubs.
}
\description{
This function calculates the proportion of shared node hubs that are connectors (HC), i.e. the top x\% of the nodes in the shared set of nodes with the highest degree that are connector nodes.
}
\details{
In this package, a tripartite network contains three groups of nodes (a-nodes,b-nodes,c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.

The function counts the proportion of connector nodes in connector node hubs (HC). The connector node hubs are the top x\% of shared nodes with the highest degree (Dominguez-Garcia and Kefi 2024). The default x\% is 20\%. It always equals 1 if all b-nodes are connector nodes.

Two types of inputs \code{network.or.subnet_mat1} can be processed:
\itemize{
\item An "igraph" object with node attribute 'level' (0 for a-nodes, 1 for b-nodes, 2 for c-nodes). If the input is a weighted network, the edge should have a 'weight' attribute.
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
hc(Net)

#empirical network
data(PPH_Coltparkmeadow)
Net <- PPH_Coltparkmeadow
hc(Net)
set.seed(13)
library(igraph)
E(Net)$weight<-runif(length(E(Net)),0.1,1)#random weights assigned
pc(Net,weighted=T)

##input as binary matrices,with row names.
set.seed(12)
md1 <- matrix(sample(c(0,1),5*8,replace=TRUE),5,8,dimnames = list(paste0("b",1:5),paste0("c",1:8)))
md2 <- matrix(sample(c(0,1),20*30,replace=TRUE),20,30,dimnames = list(paste0("b",1:20),paste0("a",1:30)))
hc(md1,md2)

##input as weighted matrices,with row numbers as row names.
set.seed(17)
mdw1 <- matrix(sample(c(rep(0,20),runif(20,0,1))),5,8)
mdw2 <- matrix(sample(c(rep(0,500),runif(100,0,1))),20,30)
hc(mdw1,mdw2)
hc(mdw1,mdw2,weighted=T)

}
\references{
Battiston, F., Nicosia, V. & Latora, V. (2014) Structural measures for
multiplex networks. Physical Review E, 89, 032804.

Domínguez-García, V. and Kéfi, S. (2024). The structure and robustness of
ecological networks with two interaction types. PLOS Computational Biology,
20(1), e1011770.

Guimera, R. & Amaral, L.A.N. (2005) Cartography of complex networks: modules
and universal roles. Journal of Statistical Mechanics: Theory and
Experiment, 2005, P02001.
}
