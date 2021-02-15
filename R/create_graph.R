create_graph <- function(spots, ids){
  ids <- ids[spots,]
  
  # create adjacency matrix
  adj.mat <- matrix(0, nrow = length(spots), ncol = length(spots))
  rownames(adj.mat) <- spots
  colnames(adj.mat) <- spots
  
  for(s in spots){
    s.x <- ids[s, "X"]
    s.y <- ids[s, "Y"]
    
    for(y in c(s.y+1, s.y-1)){
      neighbor <- intersect(which(ids[,"Y"] == y), which(ids[,"X"] == s.x))
      if(length(neighbor) > 0){
        adj.mat[s, rownames(ids)[neighbor]] <- 1
        adj.mat[rownames(ids)[neighbor], s] <- 1
      }
    }
    for(x in c(s.x+1, s.x-1)){
      neighbor <- intersect(which(ids[,"Y"] == s.y), which(ids[,"X"] == x))
      if(length(neighbor) > 0){
        adj.mat[s, rownames(ids)[neighbor]] <- 1
        adj.mat[rownames(ids)[neighbor], s] <- 1
      }
    }
  }
  
  # create graph from matrix
  neighbor.graph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.mat, 
                                                        mode = "undirected",
                                                        weighted = NULL,
                                                        add.colnames = NA,
                                                        add.rownames = NA
                                                        )
  return(neighbor.graph)
}