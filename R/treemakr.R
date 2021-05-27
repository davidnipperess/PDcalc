#' Convert a table of taxonomic relationships into a phylogenetic tree or matrix
#'
#' Converts a table of taxa into a phylogenetic tree or matrix. The phylogenetic
#' tree produced will be rooted with a topology corresponding to the
#' child-parent relationships encoded in the table and branch lengths determined
#' by the ages of taxa provided in the table. The phylogenetic matrix produced
#' records the distribution of child taxa amongst parent taxa.
#'
#' @param x is a \code{data.frame} with a row for each taxon. To produce a
#'   phylogenetic tree, every tip and internal node (including the root) must be
#'   represented by a row. \code{x} includes the following columns.
#'   \itemize{
#'   \item{\code{Taxon} is a \code{character} vector giving the names of taxa (assigned to tips and
#'   internal nodes).}
#'   \item{\code{Parent} is a \code{character} vector giving the
#'   corresponding parent taxon. Each parent name must match exactly to a name
#'   in \code{Taxon}, except for the root node which must be labelled as "Root".}
#'   \item{\code{Node_type} is a \code{character} vector labelling each row as one of
#'   "Root", "Tip" or "Internal".}
#'   \item{\code{Node_age} is \code{numeric} vector
#'   giving the age of the internal node or tip. Tips would normally have an age
#'   of 0 unless extinct. This column is not required if generating a
#'   phylogenetic matrix only.}
#'   }
#' @param output is a \code{character} vector of length 1 indicating whether the
#'   output should be a matrix (\code{output="matrix"}) or a tree
#'   (\code{output="tree"}).
#' @details \code{treemakr} takes a table of child-parent relationships among
#'   taxa and converts to a phylogenetic tree or matrix. The table could be
#'   encoded from a published phylogenetic tree, for which no newick or nexus
#'   file is available, or from a hierarchical taxonomy. If a tree is desired,
#'   node ages are required but these could be arbitrary (e.g. representing
#'   taxonomic level).
#' @return If \code{output="matrix"}, a binary \code{matrix} is returned with a
#'   row for each tip and a column for each internal and terminal (tip) node. If
#'   a tip is descendant from that node, the matrix records that relationship as
#'   a 1, and 0 otherwise. If \code{output="tree"}, a \code{phylo} object (see
#'   \code{ape} package) is returned. If the taxa in \code{x} are not in
#'   cladewise order, the resulting tree will not plot correctly in \code{ape}
#'   (the plot function requires trees to be organised in cladewise order). This
#'   can be easily fixed by writing the tree to a Newick or Nexus file (using
#'   the \code{write.tree} or \code{write.nexus} in the \code{ape} package) and
#'   then reading the file back into R.
#' @export
#' @examples 
#' data(bandicoot_table)
#' plot(treemakr(bandicoot_table,output="tree"))


treemakr <- function(x, output) {
  
  tips <- which(x$Node_type=="Tip")
  root <- which(x$Node_type=="Root")
  internals <- which(x$Node_type=="Internal")
  
  if(output=="matrix") {
    # make an MRP matrix (tips as rows, nodes (inc. tips) as columns)
    
    phylomatrix <- matrix(0, nrow=length(tips), ncol=length(internals)+length(tips)+1)
    j = 0
    
    for(i in tips) {
      ancestor <- x$Parent[i]
      node <- which(x$Taxon==ancestor)
      lineage <- c(i,node)
      j = j+1
      while(node!=root){
        ancestor <- x$Parent[node]
        node <- which(x$Taxon==ancestor)
        lineage <- c(lineage,node)}
      phylomatrix[j,lineage]=1
    }
    
    rownames(phylomatrix) <- x$Taxon[tips]
    colnames(phylomatrix) <- x$Taxon
    return(phylomatrix)
  }
  
  if(output=="tree") {
    # make a phylo object
    
    n <- length(tips)
    m <- length(internals)+1
    Nnode <- m
    edge <- matrix(nrow=n+m-1,ncol=2)
    tip.label <- x$Taxon[tips]
    node.label <- c(x$Taxon[root],x$Taxon[internals])
    taxon.id <- vector(mode="numeric",length=m+n)
    taxon.id[tips] <- 1:n
    taxon.id[root] <- n+1
    taxon.id[internals] <- 2:m+n
    edge[,2] <- taxon.id[-root]
    edge[,1] <- taxon.id[match(x$Parent[-root],x$Taxon)]
    edge.length <- x$Node_age[match(edge[,1],taxon.id)]-x$Node_age[match(edge[,2],taxon.id)]
    tree <- list(edge=edge,tip.label=tip.label,Nnode=Nnode,edge.length=edge.length,node.label=node.label)
    class(tree) <- "phylo"
    return(tree)
  }
}
