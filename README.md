
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ILSM: A package for analyzing the interconnection structure of tripartite interaction networks

<!-- badges: start -->
<!--[![CRAN_Status](http://www.r-pkg.org/badges/version/ILSM)](https://cran.r-project.org/package=ILSM) -->

[![R-CMD-check](https://github.com/WeichengSun/ILSM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WeichengSun/ILSM/actions/workflows/R-CMD-check.yaml)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/WeichengSun/ILSM/graph/badge.svg)](https://app.codecov.io/gh/WeichengSun/ILSM)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-Stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html##Stable)

<!-- badges: end -->

ILSM is designed for analyzing interconnection structures including
interconnection patterns, centrality and motifs in tripartite
interaction networks.

## Installation

You can install the development version of ILSM from GitHub:

``` r
devtools::install_github("WeichengSun/ILSM")
```

## Definition of a tripartite network

For clarification in this following context, we refer to a tripartite
network as a two-subnetwork interaction network: it is composed of three
sets of nodes (a-nodes, b-nodes, and c-nodes) and two subnetworks (the P
subnetwork with links between a-nodes and b-nodes, and the Q subnetwork
with links between b-nodes and c-nodes); the b-nodes are the shared set
of nodes; connector nodes are the common nodes of both subnetworks in
the b-nodes (Fig. 1a). No intra-guild interactions are considered if not
specified.

We provide three examples to showcase the functionality of the ‘ILSM’
package and the ecological application of interconnection motif
analysis. First, we provide a worked example to show how to calculate
interconnection patterns, motifs and centralities. Second, we compare
the differences of interconnection patterns, centralities and motifs in
characterizing variation of empirical plant-herbivore-parasitoid (PHP)
and pollinator-plant-herbivore (PPH) networks. Finally, we compare the
profiles of interconnection motifs to see if they can differentiate PHP
and PPH networks.

## A worked example of analyzing interconnection structures

As a worked example, we use a published pollinator-plant-herbivore (PPH)
binary tripartite network (Villa-Galaviz et al. 2021). This PPH network
is a subset of a large hybrid network including plants, flower visitors,
leaf minors, and parasitoids from a long-term nutrient manipulation
experiment (Colt Park Meadows) located at 300 m elevation in the
Ingleborough National Nature Reserve in the Yorkshire Dales, northern
England (54°12′N, 2°21′W). It contained pollinators, plants and
herbivores, corresponding with mutualistic interactions in the
pollinator-plant subnetwork and antagonistic interactions in the
plant-herbivore subnetwork. Plants are the shared set of species between
two subnetworks. We also generate random weights to show the analysis of
interconnection structures for weighted or quntitative tripartite
networks.

### Interconnection pattern

Interconnection pattern refers a macro-scale property on how connector
nodes (species) interconnect two subnetworks. Five interconnection
patterns are supported including proportion of connector nodes (POC),
correlation of interaction degree (CoID), correlation of interaction
similarity (CoIS), participation coefficient (PC), and proportion of
connector nodes in shared node hubs (HC).

``` r
library(ILSM);library(igraph)
data(PPH_Coltparkmeadow)
E(PPH_Coltparkmeadow)$weight<-runif(length(E(PPH_Coltparkmeadow)),0.1,1)#Generating random weights for showing the weighted interconnection structures.
#proportion of connector nodes
poc(PPH_Coltparkmeadow)
#correlation of interaction degree 
coid(PPH_Coltparkmeadow)
coid(PPH_Coltparkmeadow,weighted=T)
#correlation of interaction similarity
cois(PPH_Coltparkmeadow)
cois(PPH_Coltparkmeadow,weighted=T)
#participation coefficient 
pc(PPH_Coltparkmeadow)
pc(PPH_Coltparkmeadow,weighted=T)
#proportion of connector nodes in shared node hubs
hc(PPH_Coltparkmeadow)
hc(PPH_Coltparkmeadow,weighted=T)
```

### Interconnection motif

Here, we define 48 forms of interconneciton motifs (Fig. 2). An
interconnection motif must comprise three sets of connected nodes: the
connector nodes (belonging to b-nodes), the nodes in one subnetwork
(belonging to a-nodes in the P subnetwork), and the nodes in the other
subnetwork (belonging to c-nodes in the Q subnetwork). We further
restricted each motif to contain six nodes and no intra-guild
interactions.

The 48 interconnection motifs are provided as ‘igraph’ objects by
‘Multi_motif’ in this package.

``` r
library(ILSM)
motif_names<-c("M111","M112","M113","M114","M211","M212","M213","M311","M312","M411","M121","M122-1",
       "M122-2","M122-3","M123-1","M123-2","M123-3","M123-4","M123-5","M221-1","M221-2",
       "M221-3","M222-1","M222-2","M222-3","M222-4","M222-5","M222-6","M222-7","M222-8",
       "M222-9","M321-1","M321-2","M321-3","M321-4","M321-5","M131","M132-1","M132-2",
       "M132-3","M132-4","M132-5","M231-1","M231-2","M231-3","M231-4","M231-5","M141")
mr <- par(mfrow=c(6,8),mar=c(1,1,3,1))
IM_res<-Multi_motif("all")
 for(i in 1:48){
     plot(a[[i]],
          vertex.size=30, vertex.label=NA,
          vertex.color="#D0E7ED",main=motif_names[i])
}
par(mr)
```

The 48 interconnection motifs are named named “MABC-i”: M means
“motif’,”A” is the number of a-nodes, “B” is the number of b-nodes, “C”
is the number of c-nodes and “i” is the serial number for the motifs
with the same “ABC”. The interconnection motifs are ordered by the
number of connector nodes (from 1 to 4). The numbers from 1 to 70 in
connector nodes represent the unique roles defined by motifs.

``` r
knitr::include_graphics("../man/figure/motif_ILSM.png")
```

Fig. 2. The 48 forms of interconnection motifs with 3-6 nodes and no
intra-guild interactions. Blue and grey nodes from one subnetwork, and
grey and orange nodes from the other subnetwork. Grey nodes are
connector nodes.

For a tripartite network, the ‘icmotif_count’ gives the counts of each
motif and ‘icmotif_role’ gives counts of roles.

``` r
icmotif_count(PPH_Coltparkmeadow)
icmotif_role(PPH_Coltparkmeadow)
icmotif_count(PPH_Coltparkmeadow, weighted=T)
icmotif_role(PPH_Coltparkmeadow, weighted=T)
```

### Interconnection centrality

Interconnection centrality is to measure the importance of connector
nodes in connecting two subnetworks in a tripartite network. This is
different from common centrality measures (e.g., from ‘igraph’) that
treat connectors equally with other nodes. Three metrics are provided in
this package. For binary networks, interconnection degree centrality of
each connector species is defined as the product of its degree values
from two subnetworks, interconnection betweenness centrality is defined
by the number of shortest paths from a-nodes to c-nodes going through
connector species and interconnection closeness centrality is defined by
the inverse of the sum of distances from connector species to both
a-nodes and c-nodes. For weighted networks, interaction strengths are
used in the calculation of weighted degree, shortest path, and distance.

``` r
node_icc(PPH_Coltparkmeadow)
node_icc(PPH_Coltparkmeadow,weighted=T)
```

## License

The code is released under the MIT license (see LICENSE file).

## References

Simmons, B. I., Sweering, M. J., Schillinger, M., Dicks, L. V.,
Sutherland, W. J., & Di Clemente, R. (2019). bmotif: A package for motif
analyses of bipartite networks. Methods in Ecology and Evolution, 10(5),
695-701.

Mora, B.B., Cirtwill, A.R. and Stouffer, D.B., 2018. pymfinder: a tool
for the motif analysis of binary and quantitative complex networks.
bioRxiv, 364703.

Domínguez-García, V., & Kéfi, S. (2024). The structure and robustness of
ecological networks with two interaction types. PLOS Computational
Biology, 20(1), e1011770.

Sauve, A. M., Thébault, E., Pocock, M. J., & Fontaine, C. (2016). How
plants connect pollination and herbivory networks and their contribution
to community stability. Ecology, 97(4), 908-917.

Pilosof, S., Porter, M. A., Pascual, M., & Kéfi, S. (2017). The
multilayer nature of ecological networks. Nature Ecology & Evolution,
1(4), 0101.

Domenico, M. D. 2022. Multilayer Networks: Analysis and Visualization.
Introduction to muxViz with R. . Springer, Cham.

## Citation

Manuscript is being prepared for submission and citations are currently
available at CRAN
[10.32614/CRAN.package.ILSM](https://cran.r-project.org/web/packages/ILSM/index.html)
