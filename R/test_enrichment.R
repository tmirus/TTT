#' @export
test_enrichment <- function(analysis.data, clustering, counts, split_up_down = TRUE, db_dataset = 'mmusculus_gene_ensembl', exclusive = F){
    deg.table <- analysis.data$differential_genes
    gene.table <- analysis.data$gene.cluster.table

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
        print(cl)
        # attribute genes only to the cluster with the lowest p-value
        if(exclusive){
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
        }else{
            # use for each cluster all genes with p-value below cutoff (1e-05)
            gene.list <- gene.table[, as.character(cl)]
            names(gene.list) <- rownames(gene.table)
            gene.list <- names(gene.list[which(gene.list < 1e-5)])
            if(split_up_down){
                up_down <- sapply(gene.list, function(x){
                    sign(mean(counts[clustering == cl,x]) - mean(counts[clustering != cl, x]))
                })
                enrichments[[cl]] <- list()
                up.genes <- gene.list[which(up_down > 0)]
                down.genes <- gene.list[which(up_down < 0)]
                enrichments[[cl]][["down"]] <- enrichment_test(colnames(counts), down.genes, db, go_ids)
                enrichments[[cl]][["up"]] <- enrichment_test(colnames(counts), up.genes, db, go_ids)
            }else{
                enrichments[[cl]] <- enrichment_test(colnames(counts), gene.list, db, go_ids)
            }
        }
    }

    return(enrichments)
}
