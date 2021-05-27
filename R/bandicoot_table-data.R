#' A taxonomic table of living and recently extinct bandicoots
#'
#' A table encoding child-parent relationships and ages of origination for each
#' taxon (clade) defined in a dated and fully resolved phylogenetic tree of
#' living and recently extinct bandicoots (Perameloidea). Derived from the tree
#' published in Kear, Aplin & Westerman (2016). An example of an taxonomic table
#' used as input for the \code{treemakr} function.
#'
#' @docType data
#'
#' @usage data(bandicoot_table)
#'
#' @format A \code{data.frame} with a row for each taxon and the following columns.
#'   \itemize{
#'   \item{\code{Taxon} is a \code{character} vector giving the names of taxa
#'   assigned to tips and internal nodes.}
#'   \item{\code{Parent} is a \code{character} vector giving the corresponding parent taxon, except for the root node
#'   which must be labelled as "Root".} 
#'   \item{\code{Node_type} is a \code{character} vector labelling each row as one of "Root", "Tip" or "Internal".}
#'   \item{\code{Node_age} is \code{numeric} vector giving the age (to the nearest
#'   million years) of the internal node or tip. Tips have an age of 0, even
#'   those that are recently extinct.}
#'   }
#'
#' @keywords datasets
#'
#' @references Kear BP, Aplin KP & Westerman M (2016) Bandicoot fossils and DNA
#'   elucidate lineage antiquity amongst xeric-adapted Australasian marsupials.
#'   \emph{Scientific Reports} 6: 37537.
#'   
"bandicoot_table"
