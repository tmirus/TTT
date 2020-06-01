#' create an enrichment matrix for a given ST count matrix
#' 
#' @param counts non-negative numeric matrix. rows correspond to spots, 
#' columns correspond to genes
#' @param ncores number of threads to be used for parallel execution
#' @param db_dataset ensembl dataset to use, default 'mmusculus_gene_ensembl'
#' @return list with two entries: \cr
#' 1) enrichment.list: list containing an enrichment table for each spot in the count matrix\cr
#' 2) term.matrix: numeric matrix, terms x spots containing the p-value for each term in each spot
#' @export

calculate_spot_enrichments <- function(counts, n.cores = 4, db_dataset = 'mmusculus_gene_ensembl'){
    spot.enrichments <- calculate_enrichments(counts, db_dataset, n.cores = n.cores)
    names(spot.enrichments) <- rownames(counts)

    terms <- sapply(spot.enrichments, function(x){
        weights <- as.numeric(x[, "weightFisher"])
        names(weights) <- x[, "GO.ID"]
        weights
    })
    names(terms) <- names(spot.enrichments)

    # create a vector containing the unique GO terms
    term.names <- unique(sapply(names(unlist(terms)), function(x){strsplit(x, ".", fixed = T)[[1]][2]}))

    # matrix for storing the enrichment information, terms x spots
    term.matrix <- matrix(0, nrow = length(term.names), ncol = length(spot.enrichments))
    rownames(term.matrix) <- term.names
    colnames(term.matrix) <- names(spot.enrichments)

    for(s in names(spot.enrichments)){
        temp.terms <- names(terms[[s]])
        term.matrix[temp.terms, s] <- as.numeric(terms[[s]])
    }
    # 1 means no enrichment, < 1 are the p-values
    term.matrix[term.matrix == 0] <- 1

    #cat("Dimension of the nerichment matrix: ", dim(term.matrix), "\n", sep = "")
    #cat("Fraction of unenriched terms: ", sum(term.matrix == 1) / length(as.vector(term.matrix)), "\n", sep = "")

    return(list(enrichment.list = spot.enrichments, term.matrix = term.matrix))
}