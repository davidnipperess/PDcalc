#' Algorithm to randomly resolve polytomies in a phylogenetic tree
#' 
#' Many trees used in phylogenetic ecology include polytomies (nodes with more 
#' than two descendent branches) and are the result of uncertainty about the 
#' true branching order of taxa. Polytomies result in an overestimation of 
#' Phylogenetic Diversity because only two of the members of a polytomy will 
#' have a divergence as great as the node depth.  To overcome this problem, this
#' algorithm resolves all polytomies randomly (see details) using the procedure
#' described by Rangel et al. (2015). \emph{Please note that this algorithm is
#' designed for rooted, ultrametric trees} and the resulting fully resolved tree
#' will also be ultrametric. Because trees are resolved randomly, the algorithm
#' is designed to be run multiple times to explore the range of potential
#' solutions (which could be very large).

#' @param phy is a rooted, ultrametric phylogenetic tree with branch lengths 
#'   stored as a phylo object (as in the \code{ape} package).
#' @param runs is the number of repetitions of the algorithm (to produce 
#'   multiple equally plausible trees). Default is 1.
#' @details \code{bifurcatr} takes a phylogenetic tree (rooted and ultrametric, 
#'   with branch lengths) and randomly resolves all polytomies. The algorithm 
#'   proceeds as follows: 1) randomly choose two edges from a polytomy, 2) join 
#'   the edges with a new node, 3) generate a new edge, 4) link polytomy node to
#'   new node via new edge, 5) generate new branch length (from a random uniform
#'   distribution) for new edge, 6) adjust lengths of descendent edges (to 
#'   preserve ultrametricity). The user should repeat the algorithm (using the 
#'   \code{runs} argument) to generate several equally plausible trees. The 
#'   algorithm is a modification of that descibed in Rangel et al. (2015).
#' @return If \code{runs = 1}, a phylo object (see \code{ape} package) is 
#'   returned. If \code{runs > 1}, then a multiPhylo object (see \code{ape} 
#'   package) is returned, consisting of one tree for each run of the algorithm.
#'   Because the algorithm simply appends new edges to the resulting phylo 
#'   object, it will not plot correctly in \code{ape} (the plot function 
#'   requires trees to be organised in cladewise order). This can be easily 
#'   fixed by writing the tree (or trees) to a Newick or Nexus file (using the
#'   \code{write.tree} or \code{write.nexus} in the \code{ape} package) and then
#'   reading the file back into R.
#' @references \itemize{ \item{Rangel TF, Colwell RK, Graves GR, Fucikova K, 
#'   Rahbek C, & Diniz-Filho JAF (2015). Phylogenetic uncertainty revisited: 
#'   Implications for ecological analyses. \emph{Evolution} 69: 1301–1312.}}
#' @export
#' 


bifurcatr <- function(phy,runs=1) {
  trees <- vector("list",length=runs)
  for(i in 1:runs) {
    tree <- phy
    resolves <- Ntip(tree)-Nnode(tree)-1
    for(j in 1:resolves) {
      descendent_counts <- rle(sort(tree$edge[,1]))
      polytomies <- descendent_counts$values[which(descendent_counts$lengths>2)] # these are parent nodes with more than 2 child nodes
      if(length(polytomies)>1) target_polytomy <- sample(polytomies,size=1) else target_polytomy <- polytomies
      polytomy_edges <- which(tree$edge[,1]==target_polytomy)
      target_edges <- sample(polytomy_edges,size=2)
      new_node <- max(tree$edge)+1
      tree$edge[target_edges,1] <- new_node
      new_edge <- c(target_polytomy,new_node)
      tree$edge <- rbind(tree$edge,new_edge)
      new_length <- runif(n=1,min=0,max= min(tree$edge.length[target_edges]))
      tree$edge.length <- c(tree$edge.length,new_length)
      tree$edge.length[target_edges] <- tree$edge.length[target_edges] - new_length
      tree$Nnode <- tree$Nnode+1
    }
    trees[[i]] <- tree
  }
  if(runs==1) {
    trees <- trees[[1]]
    class(trees) <- "phylo"
  }
  else {
    class(trees) <- "multiPhylo"
  }	
  return(trees)
}
