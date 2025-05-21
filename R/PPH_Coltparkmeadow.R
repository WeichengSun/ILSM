#' An empirical tripartite pollinator-plant-herbivore network
#'
#' A pollinator-plant-herbivore tripartite network from Villa-Galaviz et. al. 2020. Journal of Animal Ecology
#'
#'
#' @details
#'
#' An 'igraph' object. This is an empirical tripartite network provided in Villa-Galaviz et. al. 2020. This network has three guilds of species: pollinators, plants and herbivores. Pollinators and plants form the subnetwork with mutualistic interactions, and plants and herbivores form the subnetwork with antagonistic interactions. No intra-guild interactions. Plants are the shared set of species.
#'
#' @name PPH_Coltparkmeadow
#'
#' @usage data(PPH_Coltparkmeadow)
#' @keywords datasets
#'
#' @srrstats {G5.1} This standard belongs here.
#'
#' @references
#'
#' Villa‐Galaviz, E., S. M. Smart, E. L. Clare, S. E. Ward, and J. Memmott. 2021. Differential effects of fertilisers on pollination and parasitoid interaction networks. Journal of animal ecology 90:404-414.
#'
#' @examples
#' data(PPH_Coltparkmeadow)
#'

load("data/PPH_Coltparkmeadow.rda")
