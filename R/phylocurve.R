#' Generate a rarefaction curve of Phylogenetic Diversity
#' 
#' Calculates a rarefaction curve giving expected phylogenetic diversity (mean 
#' and variance) for multiple values of sampling effort. Sampling effort can be 
#' defined in terms of the number of individuals, sites or species. Expected 
#' phylogenetic diversity is calculated using an exact analytical formulation 
#' (Nipperess & Matsen 2013) that is both more accurate and more computationally
#' efficient than randomisation methods.
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
#' @param stepm is the size of the interval in a sequence of numbers of 
#'   individuals, sites or species to which \code{x} is to be rarefied.
#' @param subsampling indicates whether the subsampling will be by 
#'   \code{'individual'} (default), \code{'site'} or \code{'species'}. When 
#'   there are multiple sites, rarefaction by individuals or species is done by 
#'   first pooling the sites.
#' @param replace is a \code{logical} indicating whether subsampling should be done 
#'   with (\code{TRUE}) or without (\code{FALSE} - default) replacement.
#' @details \code{phylocurve} takes community data and a rooted phylogenetic 
#'   tree (with branch lengths) and calculates expected mean and variance of 
#'   Phylogenetic Diversity (PD) for every specified value of \code{m}
#'   individuals, sites or species. \code{m} will range from 1 to the total
#'   number of individuals/sites/species in increments given by \code{stepm}.
#'   Calculations are done using the exact analytical formulae (Nipperess &
#'   Matsen, 2013) generalised from the classic equation of Hurlbert (1971).
#'   When there are multiple sites in the community data and rarefaction is by
#'   individuals or species, sites are first pooled.
#' @return a \code{matrix} object of three columns giving the expected PD values
#'   (mean and variance) for each value of \code{m}
#' @importFrom ape drop.tip
#' @references \itemize{ \item{Hurlbert (1971) The nonconcept of Species 
#'   Diversity: a critique and alternative parameters. \emph{Ecology} 52: 
#'   577-586.} \item{Nipperess & Matsen (2013) The mean and variance of 
#'   phylogenetic diversity under rarefaction. \emph{Methods in Ecology & 
#'   Evolution} 4: 566-572.}}
#' @export
#' 
#' @examples
phylocurve <- function (x, phy, stepm=1, subsampling = "individual", 
                        replace = FALSE) {

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
  ### occurrence matrix. Also set up union and intercept matrices (for calculating
  ### variance)
  
  x <- as.matrix(x)
  
  if(subsampling=="species") {
    x <- colSums(x)
    x <- ifelse(x>0,1,0) # this will normally be a vector of 1's
    commphylo <- x %*% phylomatrix  # vector of species counts per branch
    intercept_mat <- t(phylomatrix) %*% phylomatrix # matrix of shared species counts per pair of branches
    union_mat <- (intercept_mat*-1) + as.vector(commphylo)
    union_mat <- t(union_mat) + as.vector(commphylo) # matrix of total species count per pair of branches
  }
  
  if(subsampling=="individual") {
    x <- colSums(x) # vector of abundances per species
    commphylo <- x %*% phylomatrix  # vector of abundances per branch
    abundphylo <- phylomatrix * x # matrix of abundances per species per branch.
    intercept_mat <- t(abundphylo) %*% phylomatrix # matrix of shared individuals counts per pair of branches
    union_mat <- (intercept_mat*-1) + as.vector(commphylo)
    union_mat <- t(union_mat) + as.vector(commphylo) # matrix of total abundance per pair of branches
  }
  
  if(subsampling=="site") {
    commphylo <- x %*% phylomatrix # matrix of species counts per branch per site
    commphylo <- ifelse (commphylo > 0, 1, 0) # matrix of incidence of branches per site
    intercept_mat <- t(commphylo) %*% commphylo # matrix of shared site counts per pair of branches
    commphylo <- colSums(commphylo) # vector of site counts per branch
    union_mat <- (intercept_mat*-1) + commphylo
    union_mat <- t(union_mat) + commphylo # matrix of site count per pair of branches
    commphylo <- t(as.matrix(commphylo)) # convert back to matrix to allow simple vectorisation
  }
  
  ### step 5: calculate rarefied PD for each value of m
  
  if(subsampling=="site") {
    N <- length(x[,1]) # this gives the no. of sites
  }
  else {
    N <- sum(x) # this gives either the no. of individuals or the no. of species
  }
  
  m <- seq(1,N,stepm)
  
  if(m[length(m)]!=N) {
    m <- c(m,N)}
  
  pdrare <- NULL
  pdvar <- NULL
  
  if (replace == TRUE) {
    for(i in m) {
      p <- 1-(1-(commphylo/N))^i
      pdrare <- c(pdrare, sum(p * phy$edge.length))
      varmat <- (1-(union_mat/N))^i
      varmat <- varmat - t(1-p) %*% (1-p)
      varmat <- varmat * phy$edge.length
      varmat <- t(varmat) * phy$edge.length
      pdvar <- c(pdvar, sum(varmat))
    }
  }
  else {
    for(i in m) {
      p <- 1-(exp(lchoose((N-commphylo),i)-lchoose(N,i)))
      pdrare <- c(pdrare, sum(p * phy$edge.length))
      varmat <- exp(lchoose((N-union_mat),i)-lchoose(N,i))
      varmat <- varmat - t(1-p) %*% (1-p)
      varmat <- varmat * phy$edge.length
      varmat <- t(varmat) * phy$edge.length
      pdvar <- c(pdvar, sum(varmat))
    }
  }
  
  ### step 6: compile output
  
  pdvar <- round(pdvar,digits=6) # this fixes a slight rounding error
  pdcurve <- cbind(m,pdrare,pdvar)
  
  return (pdcurve)
  
}