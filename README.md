# PDcalc
### An r package for Phylogenetic Diversity (PD) calculations in ecology, biogeography and conservation biology
**PDcalc** is a R package for calculating Phylogenetic Diversity (PD). Although various PD-related functions are available in other packages, *PDcalc* is intended to be a single repository for all known extensions of the PD measure - a measurement framework that can be broadly referred to as the "PD calculus" (Faith, 2013).

The package includes functions for:

  * Classical PD
  * PD rarefaction curves
  * PD resemblance measures
  * Phylogenetic endemism
  * Algorithms for conservation decision making
  * Polytomy resolving algorithm (for ultrametric trees)
  * Expected Phylogenetic Diversity
  
### Installation

To install the latest version (requires *devtools* package):

```
devtools::install_github("davidnipperess/PDcalc",build_vignettes = TRUE)
```

### Useful references
  * Asmyhr, M.G., Linke, S., Hose, G., & Nipperess, D.A. (2014). Systematic Conservation Planning for Groundwater Ecosystems Using Phylogenetic Diversity. PLoS One, 9(12), e115132. http://doi.org/10.1371/journal.pone.0115132
  * Chao, A., Chiu, C.-H., Hsieh, T.C., Davis, T., Nipperess, D. A., & Faith, D. P. (2015). Rarefaction and extrapolation of phylogenetic diversity. *Methods in Ecology and Evolution*, 6(4), 380–388. http://doi.org/10.1111/2041-210X.12247
  * Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61, 1–10.
  * Faith D.P. (2008) Threatened species and the potential loss of phylogenetic diversity: conservation scenarios based on estimated extinction probabilities and phylogenetic risk analysis. *Conservation Biology* 22(6): 1461–1470.
  * Faith, D.P. (2013). Biodiversity and evolutionary history: useful extensions of the PD phylogenetic diversity assessment framework. *Annals of the New York Academy of Sciences*, 1289(1), 69–89. http://doi.org/10.1111/nyas.12186
  * Hartmann K & Steel M (2006) Maximizing phylogenetic diversity in biodiversity conservation: Greedy solutions to the Noah’s Ark problem. *Systematic Biology* 55(4): 644-651.
  * Minh, B., Klaere, S., & Haeseler, A. (2006). Phylogenetic Diversity within Seconds. Systematic Biology, 55(5), 769–773. http://doi.org/10.1080/10635150600981604
  * Nipperess, D.A., Faith, D.P., & Barton, K. (2010). Resemblance in phylogenetic diversity among ecological assemblages. *Journal of Vegetation Science*, 21(5), 809–820. http://doi.org/10.1111/j.1654-1103.2010.01192.x
  * Nipperess, D.A., & Matsen, F.A., IV. (2013). The mean and variance of phylogenetic diversity under rarefaction. *Methods in Ecology and Evolution*, 4, 566–572. http://doi.org/10.1111/2041-210X.12042
  * Rangel, T.F., Colwell, R.K., Graves, G.R., Fučíková, K., Rahbek, C., & Diniz-Filho, J.A.F. (2015). Phylogenetic uncertainty revisited: Implications for ecological analyses. *Evolution*, 69(5), 1301–1312. http://doi.org/10.1111/evo.12644
  * Rosauer D, Laffan S.W., Crisp M.D., Donnellan S.C. & Cook L.G. 2009. Phylogenetic endemism: a new approach for identifying geographical concentrations of evolutionary history. *Molecular Ecology*, 18(19), 4061-4072. http://doi.org/10.1111/j.1365-294X.2009.04311.x
