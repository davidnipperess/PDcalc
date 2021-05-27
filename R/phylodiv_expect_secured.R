#' Expected Phylogenetic Diversity from conservation of threatened species
#'
#' Calculates expected Phylogenetic Diversity (PD) with and without successful
#' management of a list of threatened species. Note that this function uses that
#' version of PD that always includes the path to the root of the tree.
#'
#' @param phy is a rooted phylogenetic tree with branch lengths stored as a
#'   phylo object (as in the \code{ape} package).
#' @param species is an optional \code{character} vector of species names (tree
#'   will be trimmed to match). Each species must match to a tip label in
#'   \code{phy}. By default, all tip labels in \code{phy} will be included.
#' @param managed is a \code{character} vector of species names to be managed
#'   for conservation. Each name must match to a name in \code{species}. By
#'   default, all names in \code{species} will be included.
#' @param survival is a \code{numeric} vector of survival probabilities, given
#'   no management for conservation, corresponding to each species in
#'   \code{species}.
#' @param feasibility is an optional \code{numeric} vector of probabilities of
#'   successfully managing each species in \code{managed}. By default,
#'   probability of success is assumed to be 1.
#' @param success is a \code{numeric} vector of length 1, giving the probability
#'   of survival of a species if successfully managed for conservation. Value
#'   must be less than 1.
#' @details \code{phylodiv.expected.secured} takes the probability of survival
#'   of species (before and after management) and a phylogenetic tree (rooted
#'   and with branch lengths) and calculates expected Phylogenetic Diversity
#'   (PD). PD is defined as the total length of all branches spanning a set of
#'   species (Faith, 1992). Expected PD is the summed branch length expected to
#'   survive given the probability of survival of species (Witting & Loeschcke
#'   1995; Faith 2013). An optional argument (\code{feasibility}) allows for the
#'   gain of conserving any species to be down-weighted by the probability that
#'   management will be successful.
#' @return A \code{numeric} vector giving three values corresponding to:
#'   expected PD with management; expected PD without management; and
#'   the gain in expected PD with management (i.e. the difference).
#' @importFrom ape drop.tip
#' @importFrom ape nodepath
#' @references \itemize{\item{Faith DP. 1992. Conservation evaluation and
#'   phylogenetic diversity. \emph{Biological Conservation} 61: 1-10.}
#'   \item{Faith DP. 2013. Biodiversity and evolutionary history: useful
#'   extensions of the PD phylogenetic diversity assessment framework.
#'   \emph{Annals of the New York Academy of Sciences} 1289: 69-89.}
#'   \item{Witting L. & Loeschcke V. 1995. The optimization of biodiversity
#'   conservation. \emph{Biological Conservation} 71: 205-207.}}
#' @export
#' @examples 
#' data(bandicoot_tree)
#' threatened <- c("Macrotis_lagotis","Rhynchomeles_prattorum","Peroryctes_broadbenti",
#'        "Isoodon_auratus","Perameles_gunnii","Perameles_bougainville")
#' probs <- c(0.01,0.90,0.01,0.30,0.99,0.99,0.99,0.99,0.99,0.99,0.30,0.99,0.90,0.99,0.99,0.99,0.90,0.01,0.90)
#' phylodiv.expect.secured(bandicoot_tree,managed=threatened,survival=probs)
#' 

phylodiv.expect.secured <- function (phy, species=phy$tip.label, managed=species, survival, feasibility=1, success=0.95) {
  
  # trim tree to "species"
  phy <- drop.tip (phy, which(!species%in%phy$tip.label))
  
  # convert tree to a phylomatrix
  lineages <- nodepath(phy)
  lineages <- lapply(lineages,FUN= function(x) {which(phy$edge[,2]%in%x)})
  lineages <- lineages[order(phy$tip.label)]
  
  phylomatrix <- matrix(data=0,nrow=length(phy$tip.label),ncol=length(phy$edge.length))
  
  for(i in 1:length(lineages)) {
    phylomatrix[i,lineages[[i]]] <- 1
  }
  
  # combine inputs into a common sorting order
  inputs <- merge(data.frame(species,survival),data.frame(managed,feasibility),by.x="species",by.y="managed",all.x=TRUE)
  inputs$species <- as.character(inputs$species)
  
  # create management indices
  managed <- which(inputs$species%in%managed) # indices of managed species in 'species'
  
  # survival probability
  all_managed_secured <- inputs$survival
  all_managed_secured[managed] <- all_managed_secured[managed]+(success-all_managed_secured[managed])*inputs$feasibility[managed]  # realised probabilities by factoring in feasibility
  no_secured <- inputs$survival
  
  # expected PD with successful management
  probs <- phylomatrix * all_managed_secured
  probs <- 1 - apply(1 - probs, MARGIN = 2, FUN = prod)
  exp_PD_secured <- sum(probs * phy$edge.length)
  
  # expected PD without management
  probs <- phylomatrix * no_secured
  probs <- 1 - apply(1 - probs, MARGIN = 2, FUN = prod)
  exp_PD_unsecured <- sum(probs * phy$edge.length)
  
  # output
  output <- c(exp_PD_secured,exp_PD_unsecured,exp_PD_secured-exp_PD_unsecured)
  names(output) <- c("with management","without management","gain")
  return(output)
  
}
