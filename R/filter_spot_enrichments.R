#' this function filters the enrichment matrix for interesting terms based on number of spots a term is enriched in and
#' the terms specificity to one cluster
#' 
#' @param enrichment.mat numeric matrix of spot enrichments containing p-values (GO terms x spots) 
#' as produced by test_spot_enrichment()
#' @param term.names character vector, full names of GO ids that are the rownames of enrichment.mat
#' @param clustering numeric vector assigning each spot (in the same order as in matrix) to a cluster
#' @param spot.threshold numeric > 0, number of spots a GO term needs to be enriched in
#' @param specificity.theshold numeric, 0<x<1; fraction of spots that a Term is enriched in that need 
#' to be in the same cluster for the Term to be deemed interesting
#' @return list containing two enries:\cr
#' 1) enrichment.mat - filtered matrix contining p-values of enrichments test for 
#' each spot and GO term (terms x spots) \cr
#' 2) term.names - character vector containing full names of the GO ids in the enrichment matrix
#' @export

filter_spot_enrichments <- function(enrichment.mat, term.names, clustering = NULL, spot.threshold = 10, specificity.threshold = 0.75){

    # remove rare enrichments
    term.names <- term.names[-which(apply(enrichment.mat, 1, function(x){sum(x>0)<spot.threshold}))]
    enrichment.mat <- enrichment.mat[-which(apply(enrichment.mat, 1, function(x){sum(x>0)<spot.threshold})),]

    # find enrichments that are particular to one cluster
    if(!is.null(clustering)){
        clustering.vec <- c()
        for(i in 1:nrow(enrichment.mat)){
            distribution <- table(clustering[which(enrichment.mat[i,] > 0)])
            if(any(distribution / sum(distribution) > specificity.threshold)){
                clustering.vec <- c(clustering.vec, names(distribution)[which.max(distribution)])
            }else{
                clustering.vec <- c(clustering.vec, -1)
            }
        }
    }
    cat("Number of enriched terms: ", nrow(enrichment.mat), "\n", sep = "")
    cat("Number of specific terms: ", sum(clustering.vec > 0), "\n", sep="")

    enrichment.mat <- enrichment.mat[order(clustering.vec),]
    return(list(enrichment.mat = enrichment.mat, term.names = term.names, specific.terms = term.names[which(clustering.vec > 0)], specific.ids = rownames(enrichment.mat)[which(clustering.vec > 0)]))
}
