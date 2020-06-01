#' parallelized enrichment testing for ST data
#' 
#' @param counts non-negative numeric matrix. rows correspond to spots, 
#' columns correspond to genes
#' @param db_dataset ensembl dataset to use, default 'mmusculus_gene_ensembl'
#' @param ncores number of threads to be used for parallel execution
#' @return list containing an enrichment table for each spot

calculate_enrichments <- function(counts, db_dataset = 'mmusculus_gene_ensembl', n.cores = 4){
        biomartCacheClear()
        db <- useEnsembl(biomart = "ensembl", dataset = db_dataset, mirror = 'useast')
        go_ids <- getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'),
                    filters = 'external_gene_name',
                    values = colnames(counts),
                    mart = db,
                    verbose = F)
        
        cl <- makeCluster(n.cores)
        registerDoParallel(cl)
        geneverse <- colnames(counts)
        enrichments <- foreach(i = 1:nrow(counts)) %dopar% {
                message(as.character(i))
                spot <- rownames(counts)[i]
                expressed.genes <- colnames(counts)[which(counts[spot,] > 0)]
                try(enrichment_test(geneverse, expressed.genes, spot, db, go_ids))
        }
        stopCluster(cl)
        return(enrichments)
}