#' filter term matrix based on specificity and coverage
#' 
#' @param term.mat matrix containing enrichment information, terms x spots. Term not enriched in a spot have entry 0 in that column.
#' @param spot.cluster vector containing clustering information for the spots in term.mat
#' @param exclusive.threshold numeric >0, <1; if more than this fraction of spots in which a term is enriched are in the same cluster,
#' the term is kept
#' @param cluster.marker numeric >0, <1; if more than this fraction of the spots of any cluster are enriched for a certain term, 
#' that term is kept
#' @return term.mat, filtered for interesting terms
#' @export

reduce_term_mat <- function(term.mat, spot.cluster, exclusive.threshold = 0.75, cluster.marker = 0.85) {
	marker <- sapply(rownames(term.mat), function(x){
		for(cluster in unique(spot.cluster)){
			if(sum(term.mat[x, which(spot.cluster == cluster)] > 0) > cluster.marker * length(which(spot.cluster == cluster))){
				return(TRUE)
			}
			if(sum(term.mat[x, which(spot.cluster == cluster)] > 0) / sum(term.mat[x, ] > 0) > exclusive.threshold){
				return(TRUE)
			}
		}
		return(FALSE)
	})
	n.terms <- nrow(term.mat)
	term.mat <- term.mat[marker, ]
	#cat("Fraction of terms removed: ", 1-nrow(term.mat)/n.terms, "\n", sep = "")
	return(term.mat)
}