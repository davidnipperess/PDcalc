#' Pruning algorithm to maximise gains in expected Phylogenetic Diversity from
#' successfully managing threatened species
#'
#' A pruning algorithm that seeks to maximise gains in expected Phylogenetic
#' Diversity (PD) from successfully managing threatened species from a candidate
#' list. Starting with all candidate species being assumed to be successfully
#' managed, the algorithm proceeds by stepwise removal of species from the
#' candidate list that providing the smallest gain (from conservation) at that
#' step. The pruned candidate list at each step will be optimal for conserving
#' Phylogenetic Diversity for that number of species. The reverse order in which
#' species are pruned gives the priority for conservation.
#'
#' @param phy is a rooted phylogenetic tree with branch lengths stored as a
#'   phylo object (as in the \code{ape} package).
#' @param species is an optional \code{character} vector of species names (tree
#'   will be trimmed to match). Each species must match to a tip label in
#'   \code{phy}.By default, all tip labels in \code{phy} will be included.
#' @param managed is a \code{character} vector of species names to be included
#'   on the candidate list for conservation. Each name must match to a name in
#'   \code{species}. By default, all names in \code{species} will be included as
#'   candidates.
#' @param retain is an optional \code{character} vector of species names that
#'   are never pruned from the candidate list (i.e. always assumed to be
#'   conserved). Each name must match to to name in \code{managed}.
#' @param survival is a \code{numeric} vector of survival probabilities, given
#'   no management for conservation, corresponding to each species in
#'   \code{managed}.
#' @param feasibility is an optional \code{numeric} vector of probabilities of
#'   successfully managing each species in \code{managed}. By default,
#'   probability of success is assumed to be 1.
#' @param cost is an optional \code{numeric} vector of costs of conservation
#'   corresponding to each species in \code{managed}. Default value is 1 for all
#'   candidate species.
#' @param success is a \code{numeric} vector of length 1, giving the probability
#'   of survival of a species if successfully managed for conservation. Value
#'   must be less than 1.
#' @param random is a \code{logical} vector of length 1, indicating whether
#'   species should be pruned from the candidate list at random (\code{TRUE}),
#'   rather than optimally (\code{FALSE}).
#' @details \code{phyloprunr.prob} takes a phylogenetic tree (rooted and with
#'   branch lengths) and a list of candidate species for conservation and
#'   determines the optimal set of species that will maximise expected
#'   Phylogenetic Diversity (PD) for that set size, if those species were
#'   successfully managed for conservation. The algorithm determines step-wise
#'   which candidate species provides the smallest gain to expected PD and then
#'   removes it from the candidate list, adjusting the probabilities of survival
#'   of each branch in the tree as it goes. The algorithm will continue until
#'   the list of candidates is exhausted. At each step, ties are resolved by
#'   choosing at random from the available options. For this reason, the
#'   algorithm may need to be repeated to generate several equally optimal sets.
#'   An optional argument (\code{feasibility}) allows for the gain of conserving
#'   any species to be down-weighted by the probability that management will be
#'   successful. A further optional argument (\code{cost}) allows for gains to
#'   be divided by their respective costs.
#' @return A dataframe giving the names of each species pruned (in order of
#'   pruning), the expected Phylogenetic Diversity \emph{after} pruning that species,
#'   and the total cost of management (summed across managed species) \emph{after}
#'   pruning that species.
#' @importFrom ape drop.tip
#' @importFrom ape nodepath
#' @references \itemize{ \item{Minh B., Klaere S. & Haeseler A. (2006).
#'   Phylogenetic Diversity within Seconds. \emph{Systematic Biology} 55:
#'   769–773.}}
#' @export
#'
#' @examples

phyloprunr.prob <- function(phy,species=phy$tip.label,managed=species,retain=NULL,survival,feasibility=1,cost=1,success=0.95,random=FALSE) {
  
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
    inputs <- merge(data.frame(species,survival),data.frame(managed,feasibility,cost),by.x="species",by.y="managed",all.x=TRUE)
    inputs$species <- as.character(inputs$species)
  
  # create management indices
    managed <- which(inputs$species%in%managed) # indices of managed species in 'species'
    
  # storage vectors
    prune_order <- numeric(length=length(managed)-length(retain)) # vector of indices (from 'species') in order of pruning
    exp_PD_secured <- numeric(length=length(managed)+1-length(retain))
    all_managed_secured <- inputs$survival
    all_managed_secured[managed] <- all_managed_secured[managed]+(success-all_managed_secured[managed])*inputs$feasibility[managed]  # realised probabilities by factoring in feasibility
    
  # starting expected PD
    probs <- phylomatrix * all_managed_secured
    probs <- 1 - apply(1 - probs, MARGIN = 2, FUN = prod)
    exp_PD_secured[1] <- sum(probs * phy$edge.length)
    
  # starting budget
    budget <- numeric(length=length(managed)+1-length(retain))
    budget[1] <- sum(inputs$cost[managed]) # starting budget before pruning (all candidate and retained species are managed)
    
  # remove retained species from candidate drop list before running algorithm
    retain <- which(inputs$species%in%retain) # indices of retained species in 'species'
    managed <- setdiff(managed,retain)
    
  # pruning algorithm
    
    for(i in 1:(length(managed))) {
      
      losses <- numeric(length=length(managed))
      
      for(j in 1:length(managed)) {
        test <- probs
        test_lineage <- lineages[[managed[j]]]
        test[test_lineage] <- 1-((1-probs[test_lineage])*(1-inputs$survival[managed[j]])/(1-all_managed_secured[managed[j]])) # adjusted probability if species j not secured
        losses[j] <- sum((probs[test_lineage] - test[test_lineage])*phy$edge.length[test_lineage]) / inputs$cost[managed[j]]
      }
      if(random==FALSE) {
        drop <- which(losses==min(losses))
        if(length(drop)>1) { drop <- sample(drop,size=1) }
      }
      else {
        drop <- sample(1:length(losses),size=1)
      }
      prune_order[i] <- managed[drop]
      drop_lineage <- lineages[[managed[drop]]]
      probs[drop_lineage] <- 1-((1-probs[drop_lineage])*(1-inputs$survival[managed[drop]])/(1-all_managed_secured[managed[drop]]))
      exp_PD_secured[i+1] <- exp_PD_secured[i] - losses[drop]
      budget[i+1] <- budget[i] - inputs$cost[managed[drop]]
      managed <- managed[-drop]
    }
    
  # prepare output
    prune_order <- inputs$species[prune_order]
    outputs <- data.frame(c("None",as.character(prune_order)),exp_PD_secured,budget)
    colnames(outputs) <- c("pruning.order","expected.PD.secured","budget")
    return(outputs) # table giving each species pruned (in order of pruning), the expected PD AFTER pruning that species, and the total cost of management AFTER pruning that species.
}
