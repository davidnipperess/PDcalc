#' Pairwise resemblance in Phylogenetic Diversity of ecological samples
#' 
#' Calculates the pairwise resemblance (similarity or dissimilarity) in the 
#' Phylogenetic Diversity of ecological samples. A variety of PD versions of 
#' classic ecological resemblance measures are available as options. Optional 
#' abundance weighting is implemented using the method of Nipperess et al. 
#' (2010).
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
#' @param incidence is a \code{logical} indicating whether the data are to be treated 
#'   as incidence (binary presence-absence) or abundance.
#' @param method indicates the particular form of the resemblance index you wish
#'   to use. Current options are: \code{'sorensen'} (default - 2a/a+b+c), 
#'   \code{'jaccard'} (a/a+b+c), \code{'simpson'} (a/a+min(b,c)) and 
#'   \code{'faith'} (a+0.5d/a+b+c+d).
#' @param dissim is a \code{logical} indicating whether the pairwise resemblance values
#'   should be similarity or dissimilarity (default)
#' @details Takes a community data table and a rooted phylogenetic tree (with 
#'   branch lengths) and calculates the resemblance in Phylogenetic Diversity 
#'   (PD-resemblance) of all pairwise combinations of samples/sites. The 
#'   principles for calculating PD-resemblance on incidence data are discussed 
#'   by Ferrier et al. (2007). This approach has been extended to include 
#'   abundance data (Nipperess et al. 2010).
#' @return A \code{distance} object giving the PD-resemblance of all pairwise 
#'   combinations of sample/sites in \code{x}.
#' @importFrom ape drop.tip
#' @references \itemize { \item{Ferrier, Manion, Elith & Richardson. 2007. Using
#'   generalized dissimilarity modelling to analyse and predict patterns of beta
#'   diversity in regional biodiversity assessment. \emph{Diversity & 
#'   Distributions} 13: 252-264.}
#'   \item{Nipperess, Faith & Barton. 2010. 
#'   Resemblance in Phylogenetic Diversity among ecological assemblages.
#'   \emph{Journal of Vegetation Science} 21: 809-820.}}
#' @export
#' 
#' @examples
phyloresembl <- function (x, phy, incidence = TRUE, method = "sorensen", dissim = TRUE) {

METHODS <- c("sorensen", "jaccard", "simpson", "faith")
method <- pmatch(method, METHODS)

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

### step 4: creating a community phylogeny matrix from a MRP matrix and an occurrence matrix

x <- as.matrix (x)
commphylo <- x %*% phylomatrix

# the above code performs simple matrix multiplication to produce a composite branch by sample matrix (the community phylogeny matrix) where each cell now contains either OTU richness or total abundance for each branch in each sample.

### step 5 (optional): converting a community phylogeny matrix into an incidence (0/1) form

if (incidence == TRUE) {
	commphylo[commphylo>1] <- 1 }

### step 6: calculating shared diversity (a)

a <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
for (i in 1:length(x[ ,1])) {
		for (j in 1:length(x[ ,1])) {
			a[i,j] = sum (phy$edge.length * pmin (commphylo[i, ], commphylo[j, ])) }
		}

### step 7: calculating unshared diversity (b)

b <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
for (i in 1:length(x[ ,1])) {
		for (j in 1:length(x[ ,1])) {
			b[i,j] = sum (phy$edge.length * (pmax(commphylo[i, ], commphylo[j, ]) - commphylo[j, ])) }
		}

### step 8: calculating unshared diversity (c)

c <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
for (i in 1:length(x[ ,1])) {
		for (j in 1:length(x[ ,1])) {
			c[i,j] = sum (phy$edge.length * (pmax(commphylo[i, ], commphylo[j, ]) - commphylo[i, ])) }
		}

### step 9: calculating absent diversity (d)

d <- matrix (data = NA, nrow = length(x[ ,1]), ncol = length(x[ ,1]))
maxabundance <- vector (mode = "numeric", length = length(commphylo[1, ]))
for (i in 1:length(commphylo[1, ])) {
	maxabundance[i] = max(commphylo[,i])
	}
for (i in 1:length(x[ ,1])) {
		for (j in 1:length(x[ ,1])) {
			d[i,j] = sum (phy$edge.length * (maxabundance - pmax(commphylo[i, ], commphylo[j,]))) }
		}

### step 10: calculating the pairwise similarity Index of samples

if (method == 1) {
        pdmat <- (2 * a)/((2 * a) + b + c)
    }
    else if (method == 2) {
        pdmat <- a/(a + b + c)
    }
    else if (method == 3) {
        pdmat <- a/(a + pmin(b,c))
    }
    else if (method == 4) {
        pdmat <- (a + (0.5 * d))/(a + b + c + d)
    }

rownames (pdmat) <- rownames (x)
colnames (pdmat) <- rownames (x)
pdmat <- as.dist (pdmat)

if (dissim == TRUE) {
	pdmat <- 1-pdmat }

return (pdmat)

}