#' this function is a wrapper for sctransform::vst 
#' 
#' @param counts - non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @return normalized count matrix in the same format as parameter counts
#' @export

normalize_counts <- function(counts){
  # normalize using sctransform
  # make sure no output is produced
  return(t(sctransform::vst(t(counts))[[1]]))
}