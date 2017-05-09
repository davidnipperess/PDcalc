#' Heuristic algorithm optimising the accumulation of Phylogenetic Diversity
#' from aggregating sites
#' 
#' A heuristic algorithm that seeks to prioritise sites for conservation by 
#' optimising the accumulation of Phylogenetic Diversity (PD). PD will acumulate
#' most rapidly when the samples are pooled in a particular sequence such that 
#' the gain in PD of each additional sample is maximised.
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
#' @param iterations is the number of iterations of the algorithm (to resolve 
#'   ties). Default is 1.
#' @details \code{phylorebelo} takes community data and a phylogenetic tree 
#'   (rooted and with branch lengths) and calculates the optimal sequence of 
#'   samples for accumulating Phylogenetic Diversity (PD). The algorithm starts 
#'   by first selecting the sample with the highest PD. It will then choose the 
#'   sample most complementary to the first (giving the largest gain in PD). 
#'   Additional samples are chosen based on their complementarity to the set 
#'   already chosen until all sites have been selected. At each step, ties are 
#'   resolved by choosing at random from the available equally complementary 
#'   options. The algorithm is an adaptation of that described by Rebelo & 
#'   Siegfried (1992). The phylogenetic version implemented here is that used by
#'   Asmyhr et al. (2014).
#' @return A list of two numeric matrices. The first matrix gives the optimised 
#'   sequence (rank) of the samples as rows, with a column for each iteration of
#'   the algorithm. The second matrix is gives the corresponding gains in PD for
#'   each sample (row) and iteration.
#' @importFrom ape drop.tip
#' @references \itemize{ \item{Asmyhr MG, Linke S, Hose G, & Nipperess DA. 2014.
#'   Systematic Conservation Planning for Groundwater Ecosystems Using 
#'   Phylogenetic Diversity. \emph{PLoS One} 9: e115132} \item{Rebelo AG & 
#'   Siegfried WR. 1992. Where should nature reserves be located in the Cape 
#'   Floristic Region, South Africa? Models for the spatial configuration of a 
#'   reserve network aimed at maximizing the protection of floral diversity. 
#'   \emph{Conservation Biology} 6: 243–252.}}
#' @export
#' 
#' @examples

phylorebelo <- function (x, phy, iterations=1) {

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

  ### step 6: run the heuristic
  
  sequence_matrix <- matrix(nrow=nrow(x),ncol=iterations)
  rownames(sequence_matrix) <- rownames(x)
  gains_matrix <- matrix(nrow=nrow(x),ncol=iterations)
  rownames(gains_matrix) <- rownames(x)

  for(i in 1:iterations) {
    PD <- phylomatrix %*% phy$edge.length
    rank_PD <- rank(PD,ties.method="random") # because of this randomness, need to repeat a large no. of times
    winner <- which(rank_PD==max(rank_PD))
    
    gains <- PD[winner]
    sequence <- winner
    
    for(j in 1:(nrow(x)-1)) {
      reserved <- phylomatrix[sequence,]
      if(length(sequence)>1) {
        reserved <-colSums(reserved)
      }
      complement <- phylomatrix
      complement[,which(reserved>0)] <- 0
      PD <- complement %*% phy$edge.length
      PD[sequence] <- -1
      rank_PD <- rank(PD,ties.method="random")
      winner <- which(rank_PD==max(rank_PD))
      sequence <- c(sequence,winner)
      gains <- c(gains,PD[winner])
    }
    
    sequence_matrix[sequence,i] <- 1:nrow(x)
    gains_matrix[sequence,i] <- gains
  }
  
  ### step 7: output the matrices

  return(list(sequence_matrix,gains_matrix))
}
