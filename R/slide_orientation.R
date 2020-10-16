#' this function changes the orientation of the data by
#' mirroring x or y coordinates and/or rotating to the right
#' by multiples of 90 degrees
#' 
#' @param ids barcode data frame as returned by process_input
#' @param flip.x logical, should x coordinates be flipped? default FALSE
#' @param flip.y logical, should y coordinates be flipped? default FALSE
#' @param rotate integer, multiple of 90; degrees by which data is rotated to the right
#' @example 
#' # flip slide upside down
#' ids <- slide_orientation(ids, flip.x = FALSE, flip.y = TRUE, rotate = 0)
#' @details If more than one action is specified, the tasks are carried out in the following order:\cr
#' - flip x\cr
#' - flip y\cr
#' - rotate\cr
#' @return barcode data frame
#' @export

slide_orientation <- function(ids, flip.x = FALSE, flip.y = FALSE, rotate = 0){
  # x flip
  if(flip.x){
    ids$X <- (max(ids$X) + 1) - ids$X 
  }
  # y flip
  if(flip.y){
    ids$Y <- (max(ids$Y) + 1) - ids$Y
  }
  # rotate
  if(rotate %% 90){
    for(i in 1:(rotate/90)){
      temp <- ids$Y
      ids$Y <- ids$X
      ids$X <- -temp
    }
  }else{
    warning("rotate must be a divisible by 90")
  }
  return(ids)
}