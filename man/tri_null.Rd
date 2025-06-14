% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tri_null.R
\name{tri_null}
\alias{tri_null}
\title{Generates null models for tripartite network}
\usage{
tri_null(
  trinet,
  null_N = 100,
  null_type = c("sauve", "sub_P", "sub_Q", "both_sub"),
  sub_method
)
}
\arguments{
\item{trinet}{An 'igraph' object that represents a tripartite network.}

\item{null_N}{The number of null models to be generated.  Default to 100.}

\item{null_type}{Character. Should be one of "sauve","sub_P","sub_Q" and "both_sub".See details.}

\item{sub_method}{The method to shuffle subnetworks. Must be provided when null_type ="sub_P","sub_Q" or "both_sub" .  a character specifying the null model algorithm listed on the help page of vegan::commsim. If null_type = 'sauve', it will be ignored.}
}
\value{
Return a list of null models of a tripartite network.
}
\description{
A wrapper function to generate different null models for binary and quantitative tripartite networks
}
\details{
In this package, a tripartite network contains three groups of nodes (a-nodes, b-nodes, c-nodes)  and two subnetworks (P includes the links between a-nodes and b-nodes, Q includes the links between b-nodes and c-nodes). Connector nodes belong to b-nodes.

\strong{null_type}
\itemize{
\item “sub_P”, “sub_Q” and "both_sub" use the null model algorithms from \strong{vegan::commsim} to shuffle single subnetwork or both of them.
\item “sauve” rearranges connector species without changing subnetworks, following sauve et al.(2016).
}
}
\examples{

## generate a random tripartite network
set.seed(12)
Net <- build_toy_net(11,15,16,0.2)

data(PPH_Coltparkmeadow)
Net <- PPH_Coltparkmeadow

set.seed(123)
tri_null_list<-tri_null(Net,null_type="both_sub",sub_method="r00")
set.seed(123)
tri_null_list<-tri_null(Net,null_type="sauve")
}
\references{
Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917.
}
