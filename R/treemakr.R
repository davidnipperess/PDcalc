## Function to convert a database table into useful objects in R for further analyses

  # x is a data frame with a row for each tip and internal node in a phylogenetic tree
  # columns of x must include:
    # $Taxon - name of tip or internal node
    # $Parent - name of parent node for each tip and internal node. If the root node, parent should be "Root"
    # $Node_type - coded as one of "Root", "Tip" or "Internal"
    # $Node_age - age assigned to a tip or node. Tips would normally have an age of 0 unless they are extinct.

treemakr <- function(x, output = c("matrix","tree")) {
  
  tips <- which(x$Node_type=="Tip")
  root <- which(x$Node_type=="Root")
  internals <- which(x$Node_type=="Internal")
  
  if(output=="matrix") {
    # make an MRP matrix (tips as rows, nodes (inc. tips) as columns)
    
    phylomatrix <- matrix(0, nrow=length(tips), ncol=length(internals)+length(tips)+1)
    j = 0
    
    for(i in tips) {
      ancestor <- x$Parent[i]
      node <- which(x$Taxon==ancestor)
      lineage <- c(i,node)
      j = j+1
      while(node!=root){
        ancestor <- x$Parent[node]
        node <- which(x$Taxon==ancestor)
        lineage <- c(lineage,node)}
      phylomatrix[j,lineage]=1
    }
    
    rownames(phylomatrix) <- x$Taxon[tips]
    colnames(phylomatrix) <- x$Taxon
    return(phylomatrix)
  }
  
  if(output=="tree") {
    # make a phylo object
    
    n <- length(tips)
    m <- length(internals)+1
    Nnode <- m
    edge <- matrix(nrow=n+m-1,ncol=2)
    tip.label <- x$Taxon[tips]
    node.label <- c(x$Taxon[root],x$Taxon[internals])
    taxon.id <- vector(mode="numeric",length=m+n)
    taxon.id[tips] <- 1:n
    taxon.id[root] <- n+1
    taxon.id[internals] <- 2:m+n
    edge[,2] <- taxon.id[-root]
    edge[,1] <- taxon.id[match(x$Parent[-root],x$Taxon)]
    edge.length <- x$Node_age[match(edge[,1],taxon.id)]-x$Node_age[match(edge[,2],taxon.id)]
    tree <- list(edge=edge,tip.label=tip.label,Nnode=Nnode,edge.length=edge.length,node.label=node.label)
    class(tree) <- "phylo"
    return(tree)
  }
}
