#' This function selects spots based on connectivity. Based on the assumption
#' that only spots with a certain number of connected neighbours are of interest, 
#' spots that are not sufficiently connected to the main part of the tissue are removed.
#' 
#' @param ids ids data frame as produced by process_input(), reduced to rows corresponding to the
#' spots to be retained based on image analysis (remove_background)
#' @return new ids data frame

clean_ids <- function(ids) {
  to.remove <- c()
  for(i in 1:nrow(ids)) {
    rmv <- TRUE
    x <- ids[i,2]
    y <- ids[i,1]
    neighbours <- intersect(which(ids[,2] %in% ((x-1):(x+1))), which(ids[,1] %in% ((y-1):(y+1))))
    if(length(neighbours) > 4) {
      rmv = FALSE
      next
    }else{
      for(n in neighbours){
        x <- ids[n,2]
        y <- ids[n,1]
        n_neighbours <- intersect(which(ids[,2] %in% ((x-1):(x+1))), which(ids[,1] %in% ((y-1):(y+1))))
        if(length(n_neighbours) > 5) {
          rmv = FALSE
          break
        }
      }
    }
    if(rmv) to.remove <- c(to.remove, i)
  }
  if(length(to.remove)>0) ids <- ids[-to.remove,]
  return(ids)
}