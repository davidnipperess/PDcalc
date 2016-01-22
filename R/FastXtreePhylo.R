# A fast and efficient version of the function  Xtree
# originally coded by Jens Schumacher on 17 January 2003
# and used by Petchey and Gaston (2006) to implement Petchey's
# approach to tree-based Functional Diversity.
#
# The original function by Jens Schumacher is available
# from Owen Petchey's website at:
#     http://www.thetrophiclink.org/resources/calculating-fd
# - last accessed 28 January 2014
#
# I did develop a similar version in early 2013 but
# this was unfortunately lost on 17 October 2013.
#
# Performance is improved greatly using this version
# of Schumacher's xtree algorithm. At least a halving of
# run time was found during testing and development of
# this function (usually around a 3-fold increase in speed.)
# This performance gain will be very useful when large
# community datasets are being analysed.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# See http://www.gnu.org/licenses/ for a copy of the license.
# 
# Peter D. Wilson 27-28 Jan 2014
#
# Adapted to work with objects of class "phylo" defined in the
# package ape. From the old FastXtree algorithm we borrow the principle
# that when we descend from tip (==terminal node) level towards the root,
# row entries in the currently focussed edge (ie col) of H1 are equal
# to the sum of all H1 columns representing edges above the focal edge
# in the tree. Using the rowSums() function vectorises this step making
# it very fast.
#
# PDW 10 March 2014

FastXtreePhylo <- function(thisTree)
{
  # In edge matrix in a phylo object, OTUs are indexed 1:Ntip, the root node
  # has index Ntip + 1, and internal nodes are indexed from root.node to
  # total.nodes. This is clever because it lets us very efficiently cross-reference between
  # edge indicies and the left and right terminal node indices of each edge. So,
  # compute the key values:
  root.node <- Ntip(thisTree) + 1
  total.nodes <- max(thisTree$edge)
  
  # Matrix H1 of Schumacher's Xtree method, which is exactly the same as Nipper's phylomatrix:
  H1 <- matrix(0,Ntip(thisTree),Nedge(thisTree),dimnames=list(thisTree$tip.label,NULL))
  
  # A short-cut: we can instantly populate H1 with terminal edges:
  rc.ind <- cbind("row"=1:Ntip(thisTree),"col"=which(thisTree$edge[,2] < root.node))
  H1[rc.ind] <- 1
  
  # Make a vector of internal node indices in descending order so we can traverse
  # tree from first nodes below tip nodes down to the root:
  internal.nodes <- seq(total.nodes,root.node,-1)
  
  # Now, visit each internal node accumulating subtended child edge records:
  for (thisNode in internal.nodes)
  {
    # Find the edge which has thisNode on its left:
    next.edge <- which(thisTree$edge[,2]==thisNode)
    
    # Which edges are children of the node to the right:
    child.edges <- which(thisTree$edge[,1]==thisNode)
    
    # Do the magic rowSums() trick to accumulate edges subtended by the current edge:
    H1[,next.edge] <- rowSums(H1[,child.edges])
  }
  
  # All done...package results as a named list:
  return(list("h2.prime"=thisTree$edge.length,"H1"=H1))  
}
