% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/null_model.R
\name{null_model}
\alias{null_model}
\title{Null model of tripartite network}
\usage{
null_model(
  network,
  number = NULL,
  null_type = c("subnetwork1", "subnetwork2", "all", "Savue")
)
}
\arguments{
\item{network}{A tripartite network of 'igraph' class. The network contains three groups of species and interactions within subnetwork without interactions between each group of species.}

\item{number}{A numeric value. The number of null model.  Default to NULL representing number 1.}

\item{null_type}{Logical. Four matrix scrambling algorithms. If null_type = NULL, default to "all".}
}
\value{
Return a list contains one or more elements. Each element represent a null model of tripartite network.
}
\description{
The null model could be generated according to different matrix scrambling algorithms for interconnection patterns in the tripartite network
}
\details{
\strong{null_type}

\itemize{
\item For each of the four types of null models, there are corresponding algorithms. The first type, “subnetwork1”, involved scrambling the adjacency matrix of the first and second groups of the tripartite network.
\item The second type, “subnetwork2”, focused on scrambling the adjacency matrix of the second and third groups.
\item Comprehensively, the third type, “all”, blended the approaches of the first two to disarrange the entire network's adjacency matrix, achieving a thorough perturbation of the network's structure.
\item The last type named “Savue” that disarranged inherent structure in terms of the groups of species connected by each interconnecting species of every subnetworks, thus exhibiting different interconnection patterns.
}

\strong{network}

About a network of type "igraph", It can be obtained from the connection matrices of subnetworks by the function \code{igraph_from_matrices}
}
\examples{

## generate a random tripartite network
set.seed(12)
Net <- build_net(11,15,16,0.2)

data(PPH_Coltparkmeadow)
Net <- PPH_Coltparkmeadow

set.seed(123)
null_model(Net)
set.seed(123)
null_model(Net,null_type="subnetwork1")
set.seed(123)
null_model(Net,null_type="Savue")
set.seed(123)
null_model(Net,number=2,null_type="Savue")
}
\references{
Vázquez, D. P., C. J. Melian, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction strength in ecological networks. Oikos 116: 1120-1127.

Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
}
