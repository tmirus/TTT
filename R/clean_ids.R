#' This function selects spots based on connectivity. Based on the assumption
#' that only spots with a certain number of connected neighbours are of interest, 
#' spots that are not sufficiently connected to the main part of the tissue are removed.
#' 
#' @param ids ids data frame as produced by process_input(), reduced to rows corresponding to the
#' spots to be retained based on image analysis (remove_background)
#' @return new ids data frame

clean_ids <- function(ids) {
  to.remove <- c()
  # iterate over all spots in ids
  for(i in 1:nrow(ids)) {
    rmv <- TRUE
    x <- ids[i,"X"]
    y <- ids[i,"Y"]
    # find spots connected to spot i
    neighbours <- intersect(
      which(ids[,"X"] %in% ((x-1):(x+1))), 
      which(ids[,"Y"] %in% ((y-1):(y+1)))
    )
    # three (+ self) directly connected spots are enough to keep
    # the current spot
    if(length(neighbours) > 4) {
      rmv = FALSE
      next
    }else{
      if(length(neighbours) > 0){
      # any neighbour must have 4 (+ self) direct neighbours,
      # otherwise the spot will be removed
      for(n in neighbours){
        x <- ids[n,"X"]
        y <- ids[n,"Y"]
        n_neighbours <- intersect(
          which(ids[,"X"] %in% ((x-1):(x+1))), 
          which(ids[,"Y"] %in% ((y-1):(y+1)))
        )
        if(length(n_neighbours) > 5) {
          rmv = FALSE
          break
        }
      }
      }
    }
    if(rmv) 
	    to.remove <- c(to.remove, i)
  }
  if(length(to.remove)>0) 
	  ids <- ids[-to.remove,]
  return(ids)
}
