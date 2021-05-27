#' Phylogenetic Endemism of ecological samples
#' 
#' Calculates phylogenetic endemism (sum of 'unique' branch lengths) of multiple
#' ecological samples.
#' 
#' @param x is the community data given as a \code{data.frame} or \code{matrix} 
#'   with species/OTUs as columns and samples/sites as rows (like in the 
#'   \code{vegan} package). Columns are labelled with the names of the 
#'   species/OTUs. Rows are labelled with the names of the samples/sites. Data 
#'   can be either abundance or incidence (0/1). Column labels must match tip 
#'   labels in the phylogenetic tree exactly!
#' @param phy is a rooted phylogenetic tree with branch lengths stored as a 
#'   phylo object (as in the \code{ape} package) with terminal nodes labelled 
#'   with names matching those of the community data table. Note that the 
#'   function trims away any terminal taxa not present in the community data 
#'   table, so it is not necessary to do this beforehand.
#' @param weighted is a \code{logical} indicating whether weighted endemism 
#'   (default) or strict endemism should be calculated.
#' @details Takes a community data table and a rooted phylogenetic tree (with 
#'   branch lengths) and calculates either strict or weighted endemism in 
#'   Phylogenetic Diversity (PD). Strict endemism equates to the total amount of
#'   branch length found only in the sample/site and is described by Faith et 
#'   al. (2004) as PD-endemism. Weighted endemism calculates the "spatial 
#'   uniqueness" of each branch in the tree by taking the reciprocal of its 
#'   range, multiplying by branch length and summing for all branch lengths 
#'   present at a sample/site. Range is calculated simply as the total number of
#'   samples/sites at which the branch is present. This latter approach is 
#'   described by Rosauer et al. (2009) as \emph{Phylogenetic endemism}.
#' @return A \code{vector} object giving the phylogenetic endemism of all 
#'   sample/sites in \code{x}.
#' @importFrom ape drop.tip
#' @references \itemize{ \item{Faith, Reid & Hunter. 2004. Integrating 
#'   phylogenetic diversity, complementarity, and endemism for conservation 
#'   assessment. \emph{Conservation Biology} 18(1): 255-261.} \item{Rosauer, 
#'   Laffan, Crisp, Donnellan & Cook. 2009. Phylogenetic endemism: a new
#'   approach for identifying geographical concentrations of evolutionary 
#'   history. \emph{Molecular Ecology} 18(19): 4061-4072.}}
#' @export
#' 

phyloendemism <- function (x, phy, weighted=TRUE) {

library(ape)

### step 1: trimming the tree to match the community data table thus creating a "community tree" (sensu Cam Webb).

  taxon_check <- phylomatchr(x,phy)
  
  if(length(taxon_check[[2]])>0) {
    stop('The following taxa in the community data table were not found on the',
         ' tips of the tree: ', paste(taxon_check[[2]], collapse=', '),
         '.\nPlease fix your taxonomy!',
         call.=FALSE)
  }
  
  if(length(taxon_check[[1]])>0) {
    warning(length(taxon_check[[1]]), " taxa were not in x, and were trimmed",
            " from the tree prior to calculation of PD.")
    phy <- drop.tip (phy, taxon_check[[1]])
  }

### step 2: converting a community tree into a MRP matrix

  # A MRP matrix, used in supertree methods, is where the membership of an OTU
  # in a clade spanned by a branch is indicated by a 0 (no) or 1 (yes). Unlike 
  # supertree MRP matrices, our matrix includes terminal branches. using Peter 
  # Wilson's function
  
  phylomatrix <- FastXtreePhylo(phy)$H1
  
  # this script fills a matrix of 0's with 1's corresponding to the incidence of
  # an OTU on a particular branch.

### step 3: re-ordering the OTUs to match.

phylomatrix <- phylomatrix[order(phy$tip.label), ]
x <- x[ ,order(colnames(x))]

# data are sorted to a common ordering standard, that is alphabetic order, so that OTUs match up.

### step 4: creating a community phylogeny matrix indicating incidence of branches per site

x <- as.matrix (x)
commphylo <- x %*% phylomatrix
commphylo <- ifelse (commphylo > 0, 1, 0)

### step 5: calculate the range (in sites) of each branch

ranges <- colSums(commphylo)

### step 6: calculate weightings per branch

weights <- t(commphylo) / ranges
if(weighted==F) weights <- ifelse(weights<1,0,1)

### step 6: calculating PD per sample

pd <- t(weights) %*% phy$edge.length

return (pd)

}