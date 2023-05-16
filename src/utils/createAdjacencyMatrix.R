create.adjacencyMatrix <- function(nNodes, edges) {
  edgesPerNode <- sapply(edges,length)
  cumEdgesPerNode <- c(0,edgesPerNode)
  for (i in 1:(length(cumEdgesPerNode)-2))
    cumEdgesPerNode[(i+2):length(cumEdgesPerNode)] <- cumEdgesPerNode[(i+2):length(cumEdgesPerNode)] + rep(edgesPerNode[[i]],length(cumEdgesPerNode[(i+2):length(cumEdgesPerNode)]))
  nEdges <- sum(edgesPerNode)
  edgeMatrix <- matrix(ncol=2,nrow=nEdges)
  for (j in 1:length(edges))
    if (length(edges[[j]]) != 0)
      for (k in 1:length(edges[[j]]))
        edgeMatrix[k+cumEdgesPerNode[j],] <- c(j,edges[[j]][k])
  adj <- matrix(data=0,nrow=nNodes,ncol=nNodes)
  for (i in 1:nrow(edgeMatrix))
    adj[edgeMatrix[i,1],edgeMatrix[i,2]] <- 1
  
  return(adj)
}