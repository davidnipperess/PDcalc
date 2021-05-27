#' Phylogenetic Diversity of ecological samples
#' 
#' Calculates the Phylogenetic Diversity (PD) of multiple samples 
#' simultaneously. Note that this function uses that version of PD that always 
#' includes the path to the root of the tree.
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
#' @details \code{phylodiv} takes community data and a phylogenetic tree (rooted
#'   and with branch lengths) and calculates the Phylogenetic Diversity (PD) of 
#'   all samples/sites. PD is defined as the total length of all branches 
#'   spanning a set of terminal taxa representing an ecological sample (Faith, 
#'   1992). Please note that, if the common ancestor (node) of the set of taxa 
#'   of a sample is not the root of the tree, then the set of branches 
#'   connecting this node to the root are also included in the calculation. 
#'   Calculations are achieved using the efficient matrix algebra solution of 
#'   Rodrigues & Gaston (2002).
#' @return A vector of Phylogenetic Diversity (PD) values for each sample/site 
#'   in \code{x}.
#' @importFrom ape drop.tip
#' @references \itemize{\item{Faith DP. 1992. Conservation evaluation and 
#'   phylogenetic diversity. \emph{Biological Conservation} 61: 1-10.} 
#'   \item{Rodrigues A & Gaston KJ. 2002. Maximising phylogenetic diversity in 
#'   the selection of networks of conservation areas. \emph{Biological
#'   Conservation} 105: 103-111.}}
#' @export
#' 

phylodiv <- function (x, phy) {
  
  ### step 1: matching taxa and trimming the tree to match the community data
  ### table thus creating a "community tree" (sensu Cam Webb).
  
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
  
  ### step 3: re-ordering the OTUs of the occurrence table and MRP matrix to
  ### match.
  
  phylomatrix <- phylomatrix[order(phy$tip.label), ]
  x <- x[ ,order(colnames(x))]
  
  # data are sorted to a common ordering standard, that is alphabetic order, so
  # that OTUs match up.
  
  ### step 4: creating a community phylogeny matrix from a MRP matrix and an
  ### occurrence matrix
  
  x <- as.matrix(x)
  phylomatrix <- x %*% phylomatrix
  
  # the above code performs simple matrix multiplication to produce a composite
  # branch by sample matrix (the community phylogeny matrix) where each cell now
  # contains either OTU richness or total abundance for each branch in each
  # sample.
  
  ### step 5: converting a community phylogeny matrix into an incidence (0/1)
  ### form
  
  phylomatrix[phylomatrix>1] <- 1
  
  ### step 6: calculating PD per sample
  
  pd <- phylomatrix %*% phy$edge.length
  
  return (pd)
  
}