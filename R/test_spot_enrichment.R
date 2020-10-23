#' perform enrichment testing based on differential gene expression analysis
#' 
#' @param analysis.data output of analyze_clustering
#' @param counts numeric non-negativ matrix containing gene expression for all spots (spots x genes)
#' @param db_dataset ensembl dataset to use. default "mmusculus_gene_ensembl" (mouse), for human data use "hsapiens_gene_ensembl"
#' @param gene_id character denoting the gene symbols to use; default "hgnc" for HGNC symbols; use "ensembl" for ensembl gene ids
#' @param ncores integer > 0, number of cores to use
#' @param mirror string indicating ensembl mirror site; must be one of "standard", "uswest", "useast", "asia"; default is "standard"
#' @return list containing two entries:\cr
#' 1) enrichments - matrix contining p-values of enrichments test for 
#' each spot and GO term enriched in at least one spot(terms x spots) \cr
#' 2) term.names - character vector containing full names of the GO ids in the enrichment matrix
#' @export

test_spot_enrichment <- function(analysis.data, counts, db_dataset = 'mmusculus_gene_ensembl', gene_id = "hgnc", ncores = 4, mirror = NULL){
    # avoid ssl certificate error when accessing ensembl
    suppressMessages(library(httr, quietly = TRUE))
    httr::set_config(config(ssl_verifypeer = FALSE))

    if(!is.null(mirror)){
	    if(!mirror %in% c("useast", "uswest", "asia")){
		    mirror <- "https://www.ensembl.org"
		    warning("mirror could not be identified")
	    }
    }else{
    	if(is.null(mirror)){
	    mirror <- "https://www.ensembl.org"
    	}else{
		if(mirror == "standard"){
			mirror <- "https://www.ensembl.org"
		}else{
			mirror <- paste0("https://",mirror,".ensembl.org")
		}
	}
    }

    # required for parallelization
	suppressMessages(library(parallel, quietly = TRUE))
	suppressMessages(library(doParallel, quietly = TRUE))
	suppressMessages(library(foreach, quietly = TRUE))

    deg.table <- analysis.data$differential_genes
    gene.table <- analysis.data$gene.cluster.table

    # getBM parameter
    if(gene_id == "hgnc"){
        gene_id <- "external_gene_name"
    }else if(gene_id == "ensembl"){
        gene_id <- "ensembl_gene_id"
    }

    # get GO database for testing
    suppressWarnings(suppressMessages({
    db <- biomaRt::useEnsembl('ENSEMBL_MART_ENSEMBL',dataset=db_dataset, host="https://uswest.ensembl.org", verbose = FALSE)
    go_ids <- biomaRt::getBM(
        attributes=c('go_id', gene_id, 'namespace_1003'), 
        filters = gene_id, 
        values = colnames(counts), 
        mart = db, 
        verbose = TRUE
    )
    }))
    
    # do an enrichment test for each spot based on filtered genes 
    # that are expressed in this spot (>5)
    cl <- makeCluster(ncores)
	registerDoParallel(cl)
    enrichments <- foreach(s = rownames(counts)) %dopar% {
        temp <- try(enrichment_test(
            colnames(counts), 
            intersect(rownames(deg.table), 
                      colnames(counts)[which(counts[s,] > 5)]), 
            db, go_ids))
        if(class(temp) != "try-error"){
            return(temp)
        }else{
            return(NULL)
        }
    }
    stopCluster(cl)

    # create list of all enriched terms
    enriched.terms <- unique(unlist(sapply(enrichments, function(x){if(!is.null(x)){x$`GO.ID`}})))
    term.names <- sapply(enriched.terms, function(t){GO.db::GOTERM[[t]]@Term})

    # matrix for storing enrichment information of all terms across spots
    enrichment.mat <- matrix(0, nrow = length(enriched.terms), ncol = nrow(counts))
    rownames(enrichment.mat) <- enriched.terms
    colnames(enrichment.mat) <- rownames(counts)

    # fill matrix with test scores (weightFisher)
    for(i in seq_len(nrow(counts))){
        if(is.null(enrichments[[i]])) next
        if(any(rownames(enrichment.mat) %in% enrichments[[i]]$GO.ID)){
            enrichment.mat[enrichments[[i]]$GO.ID, i] <- as.numeric(enrichments[[i]]$weightFisher)
        }
    }

    return(list(enrichments = enrichment.mat, term.names = term.names))
}
