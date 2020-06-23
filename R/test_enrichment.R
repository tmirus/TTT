#' @export
test_enrichment <- function(deg.table, counts, split_up_down = TRUE, db_dataset = 'mmusculus_gene_ensembl'){

    # get GO database for testing
    db <- useMart('ENSEMBL_MART_ENSEMBL',dataset=db_dataset, host="www.ensembl.org", verbose = TRUE)
    go_ids <- getBM(
        attributes=c('go_id', 'external_gene_name', 'namespace_1003'), 
        filters = 'external_gene_name', 
        values = colnames(counts), 
        mart = db, 
        verbose = F
    )

    # calculate enrichments for up and downregulated genes in each cluster
    enrichments <- list()
    for(cl in unique(deg.table$cluster)){
        sub.table <- deg.table[which(deg.table$cluster == cl),]
    	if(split_up_down){
        	enrichments[[cl]] <- list()
        	up.genes <- sub.table[which(sub.table$regulation > 0), "gene"]
        	down.genes <- sub.table[which(sub.table$regulation < 0), "gene"]

        	enrichments[[cl]][["down"]] <- enrichment_test(colnames(counts), down.genes, db, go_ids)
        	enrichments[[cl]][["up"]] <- enrichment_test(colnames(counts), up.genes, db, go_ids)
	}else{
		enrichments[[cl]] <- enrichment_test(colnames(counts), sub.table$gene, db, go_ids)
	}
    }

    return(enrichments)
}
