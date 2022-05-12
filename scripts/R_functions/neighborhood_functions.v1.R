anyNA <- function(x) any(is.na(x))
library(igraph)
get.neighborhood <- function(g,nodes,mode = "all",max.layer = 10,max.prop = 0.5)
{
 ######## inputs
 # g = igraph object,
 # nodes = root nodes to expand neighborhoods from,
 # mode = choose from c("all","in","out"): specify the traversing direction
 # max.layer = the maximum number of layers to traverse. 
 neigh.layers <- vector("list",length(nodes));names(neigh.layers) <- nodes;
 for (i in 1:length(nodes))
 {
  neigh.out <- lapply(1:max.layer,function(n,g,v,mm) {out <- neighborhood(g,order = n,nodes = v,mode = mm)[[1]];out <- V(g)$name[out];return(out)},g = g,v = nodes[i],mm = mode)
  
  # identify maximum layer that traversed less than 90% of all nodes. 
  mx.layer <- max(c(1,max(which(sapply(neigh.out,length) <= (max.prop*vcount(g)))))) 
  neigh.out <- neigh.out[1:mx.layer]
  # identify maximum layer that reaches the maximum number of neighborhood
  neigh.size <- sapply(neigh.out,length)
  max.num <- min(which(neigh.size == max(neigh.size)));
  neigh.out <- neigh.out[1:max.num]
  
  neigh.layers[[i]] <- neigh.out;
  names(neigh.layers[[i]]) <- paste("n.layer",1:length(neigh.out),sep = "_")
 }
 return(neigh.layers)
} 

read.geneSet <- function(geneSet.file)
{
 gene.list <- readLines(geneSet.file)
 gene.list <- strsplit(gene.list,"\t")
 
 names(gene.list) <- sapply(gene.list,function(x) x[1])
 gene.list <- lapply(gene.list,function(x) x[2:length(x)])
 return(gene.list)
}

output.geneSet.file <- function(geneSet,outputfname)
{
 if (!is.list(geneSet)) stop("geneSet is not a list.")
 if (is.null(names(geneSet))) stop("names(geneSet) is not defined properly.")
  
 sink(outputfname)
 cat(paste(paste(names(geneSet),"\t",sapply(geneSet,function(x) paste(x,collapse = "\t")),sep = ""),collapse = "\n"))
 sink()

 return(0)
}
