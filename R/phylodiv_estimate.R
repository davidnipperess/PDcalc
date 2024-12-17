#' Estimate Phylogenetic Diversity from a sample
#'
#' Estimates total phylogenetic diversity (both detected and undetected) for
#' each sample/site or a set of samples/sites. Estimated phylogenetic diversity
#' is calculated using an extrapolation of the rarefaction curve (Chao et al.
#' 2015). For undetected phylogenetic diversity, this is considered a lower
#' bound.
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
#' @param subsampling indicates whether the subsampling will be by
#'   \code{'individual'} (default) or \code{'site'}.
#' @details \code{phylodiv_estimate} takes community data and a rooted
#'   phylogenetic tree (with branch lengths) and calculates estimated total
#'   Phylogenetic Diversity (PD) for both detected and undetected species. When
#'   there are multiple sites in the community data and subsampling is by site,
#'   data are converted to incidence (if not already) and then pooled. When
#'   subsampling is by individual, estimated PD is calculated for every
#'   sample/site.
#' @return a \code{numeric} vector giving estimated PD values for every
#'   sample/site or a single value for pooled samples/sites.
#' @importFrom ape drop.tip
#' @references \itemize{ \item{Chao et al. (2015) Rarefaction and extrapolation
#'   of phylogenetic diversity. \emph{Methods in Ecology & Evolution} 6:
#'   380-388.}}
#' @export

phylodiv.estimate <- function (x, phy, subsampling = "individual") {
  
  ### step 1: checking if taxa match and, if so, trimming the tree to match the
  ### community data table thus creating a "community tree" (sensu Cam Webb).
  
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
  
  # A MRP matrix, used in supertree methods, is where the membership of an OTU in
  # a clade spanned by a branch is indicated by a 0 (no) or 1 (yes). Unlike
  # supertree MRP matrices, our matrix includes terminal branches.
  
  # using Peter Wilson's fast function
  phylomatrix <- FastXtreePhylo(phy)
  phylomatrix <- phylomatrix$H1
  
  ### step 3: re-ordering the OTUs of the occurrence table and MRP matrix to match.
  
  phylomatrix <- phylomatrix[order(phy$tip.label), ]
  x <- x[ ,order(colnames(x))]
  
  # data are sorted to a common ordering standard, that is alphabetic order, so
  # that OTUs match up.
  
  ### step 4: creating a community phylogeny matrix from a MRP matrix and an
  ### occurrence matrix.
  
  x <- as.matrix(x)
  
  if(subsampling=="individual") {
    commphylo <- x %*% phylomatrix  # matrix of abundances per branch per site
  }
  
  if(subsampling=="site") {
    x <- ifelse(x>0,1,0) # convert to incidence
    x <- colSums(x) # vector of site frequencies of each species/OTU
    commphylo <- x %*% phylomatrix # vector of site frequencies per branch
  }
  
  ### step 5: calculate estimated PD for each site (or pool of sites)
  
  if(subsampling=="site") {
    n <- length(x) # this gives the no. of sites
    f1 <- length(which(x==1)) # count of singleton species (found in only 1 site)
    f2 <- length(which(x==2)) # count of doubleton species (found in 2 sites)
    g1 <- sum(phy$edge.length[which(commphylo==1)]) # sum of singleton branches
    g2 <- sum(phy$edge.length[which(commphylo==2)]) # sum of doubleton branches
    if(g2>(g1*f2)/(2*f1)) {
      g0 <- ((n-1)/n) * g1^2/(2*g2)
    }
    else {
      g0 <- ((n-1)/n) * (g1*(f1-1))/(2*(f2+1))
    }
    PD_obs <- sum(ifelse(commphylo>0,1,0)*phy$edge.length) # observed PD from pooled sites
  }
  if(subsampling=="individual") {
    n <- rowSums(x) # no. of individuals per site
    f1 <- rowSums(ifelse(x==1,1,0)) # vector of counts of singleton species per site
    f2 <- rowSums(ifelse(x==2,1,0)) # vector of counts of doubleton species per site
    g1 <- ifelse(commphylo==1,1,0) # matrix of branches per site that are singletons
    g1 <- g1 %*% phy$edge.length # vector of summed lengths for singleton branches per site
    g1 <- g1[,1]
    g2 <- ifelse(commphylo==2,1,0) # matrix of branches per site that are doubletons
    g2 <- g2 %*% phy$edge.length # vector of summed lengths for doubleton branches per site
    g2 <- g2[,1]
    g0 <- numeric(length=length(n))
    standardChao <- which(g2>(g1*f2)/(2*f1))
    g0[standardChao] <- ((n[standardChao]-1)/n[standardChao]) * g1[standardChao]^2/(2*g2[standardChao])
    altChao <- which(!g2>(g1*f2)/(2*f1))
    g0[altChao] <- ((n[altChao]-1)/n[altChao]) * (g1[altChao]*(f1[altChao]-1))/(2*(f2[altChao]+1))
    PD_obs <- ifelse(commphylo>0,1,0) %*% phy$edge.length # observed PD per site
  }
  
  ### step 6: compile output
  
  PD_est <- PD_obs + g0
  return (PD_est)
  
}