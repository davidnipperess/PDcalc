#' Pruning algorithm to maximise Phylogenetic Diversity for a set of species of 
#' a given size
#' 
#' A pruning algorithm that seeks to maximise Phylogenetic Diversity (PD) for a 
#' set of species of a given size. The algorithm proceeds in a step-wise 
#' fashion, removing species that contribute the least to PD until the desired 
#' set size is achieved.
#' @param phy is a rooted phylogenetic tree with branch lengths stored as a 
#'   phylo object (as in the \code{ape} package).
#' @param size is the number of tips (species/OTUs) to prune the tree to.
#' @param iterations is the number of iterations of the algorithm (to resolve 
#'   ties). Default is 1.
#' @param trees is a \code{logical}, indicating whether sub-trees (\code{TRUE} -
#'   default) or only species lists/tip labels (\code{FALSE}) should be returned
#'   for each set that maximises PD.
#' @details \code{phyloprunr} takes a phylogenetic tree (rooted and with branch 
#'   lengths) and determines the optimal set of tips (species/OTUs) of a given 
#'   size that maximises Phylogenetic Diversity (PD). The algorithm determines 
#'   step-wise which tip contributes the least to PD (has the shortest terminal 
#'   branch) and then removes it, adjusting the terminal branch lengths of the 
#'   remaining tips (to accommodate the tip removal) as it goes. The algorithm 
#'   will continue until the stipulated set size is achieved. The resulting set 
#'   is the set that maximises PD for its size. At each step, ties are resolved 
#'   by choosing at random from the available options. For this reason, the user
#'   can repeat the algorithm (using the \code{iterations} argument) to generate
#'   several equally optimal sets. For a fully resolved ultrametric tree, there 
#'   will be 2^(n-k) possible solutions (where \emph{n} is the number of tips 
#'   and \emph{k} is the desired set size), so a large number of iterations of 
#'   the algorithm may be needed to explore the range of possible solutions. The
#'   algorithm is descibed in Minh et al. (2006).
#' @return If \code{trees = TRUE}, a multiphylo object (see \code{ape} package) 
#'   is returned, consisting of one tree for each iteration of the algorithm. If 
#'   \code{trees = FALSE}, a list of length \code{iterations} is returned where 
#'   each element is a vector of tip labels / species names corresponding to an
#'   optimal set.
#' @importFrom ape drop.tip
#' @references \itemize{ \item{Minh B., Klaere S. & Haeseler A. (2006). 
#'   Phylogenetic Diversity within Seconds. \emph{Systematic Biology} 55:
#'   769–773.}}
#' @export
#' 
#' @examples

phyloprunr <- function(phy, size, iterations=1, trees=TRUE) {
  max_PD_trees <- vector(mode="list",length=iterations)
  for(i in 1:iterations) {
    subphy <- phy
    while(Ntip(subphy)>size) {
      tips <- which(subphy$edge[,2]<(Ntip(subphy)+1))
      shortest_tip <- which(subphy$edge.length[tips]==min(subphy$edge.length[tips]))
      shortest_tip <- tips[shortest_tip]
      shortest_tip <- subphy$edge[shortest_tip,2]
      if(length(shortest_tip)>1) {
        shortest_tip <- sample(shortest_tip, size=1, replace=FALSE) # take a random tip to resolve ties
      }
      subphy <- drop.tip(subphy, tip=shortest_tip)
    }
    if(trees==TRUE) {
      max_PD_trees[[i]] <- subphy
    }
    else {
      max_PD_trees[[i]] <- subphy$tip.label # more efficient storage because topology won't change but labels might
    }
  }
  if(trees==TRUE) {
    class(max_PD_trees) <- "multiPhylo"
  }
  return(max_PD_trees)
}