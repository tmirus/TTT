#' perform enrichment testing based on differential gene expression analysis
#' 
#' @param analysis.data output of analyze_clustering
#' @param clustering numeric vector assigning clusters to all spots (same order as spots in counts matrix)
#' @param counts numeric non-negativ matrix containing gene expression for all spots (spots x genes)
#' @param split_up_down boolean, should enrichment testing be performed separately on the list of up and downregulated genes for each cluster? default TRUE
#' @param db_dataset ensembl dataset to use. default "mmusculus_gene_ensembl" (mouse), for human data use "hsapiens_gene_ensembl"
#' @param exclusive boolean, should each gene only be considered for enrichment analysis in the cluster where it has the lowest p-value? default FALSE
#' @param gene_id character denoting the gene symbols to use; default "hgnc" for HGNC symbols; use "ensembl" for ensembl gene ids
#' @return list containing enrichmnet information for each cluster. if split_up_down is TRUE, 
#' the list entry for each cluster is a list with two elements:\cr
#' 1) up - table with enrichment information for upregulated genes (as returned by topGO)\cr
#' 2) down - table with enrichment information for downregulated genes (as returned by topGO)\cr
#' If split_up_down is false, the list contains one enrichment table for each cluster.
#' @export

test_enrichment <- function(analysis.data, clustering, counts, split_up_down = TRUE, db_dataset = 'mmusculus_gene_ensembl', exclusive = F, gene_id = "hgnc"){
    deg.table <- analysis.data$differential_genes
    gene.table <- analysis.data$gene.cluster.table

    if(gene_id == "hgnc"){
        gene_id <- "external_gene_name"
    }else if(gene_id == "ensembl"){
        gene_id <- "ensembl_gene_id"
    }

    # get GO database for testing
    suppressMessages({
        db <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL',dataset=db_dataset, host="www.ensembl.org", verbose = TRUE)
        go_ids <- biomaRt::getBM(
            attributes=c('go_id', gene_id, 'namespace_1003'), 
            filters = gene_id, 
            values = colnames(counts), 
            mart = db, 
            verbose = F
        )
    })

    # calculate enrichments for up and downregulated genes in each cluster
    enrichments <- list()
    for(cl in unique(deg.table$cluster)){
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
