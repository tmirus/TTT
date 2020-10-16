#' this function returns a table of interesting genes based on lasso model fit, specificity and differential expression
#' 
#' @param differential.genes data frame containing information about differential expression of genes
#' @param entropies named numeric vector of entropies for all genes
#' @param lasso.data output of build_lassos function; default NULL
#' @param deg.weight weight of p-value in overall gene ranking; default 2
#' @param lls.weight weight of lasso fit in overall gene ranking; default 3
#' @param entropy.weight weight of entropy in overall gene ranking; default 1
#' @param build.lasso logical, if lasso.data is NULL, should lasso models be calculated? default FALSE
#' @param gamma numeric > 0, if build.lasso is TRUE, the gamma value passed to build_lassos; default 3
#' @param ncores numeric > 0, the number of cores available for calculation; only used for lasso building, default 4
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param ids data frame or matrix assigning spatial coordinates to the spots (see process_input() for details)
#' @param verbose logical, default TRUE
#' @return data frame containing gene name, cluster with lowest p-value, lowest p-value, 
#' up/downregulation information, lls fit (if lasso data was supplied) and entropy for each gene
#' @export

filter_genes <- function(differential.genes, entropies, lasso.data = NULL, deg.weight = 2, lls.weight = 3, entropy.weight = 1, build.lasso = FALSE, gamma = 3, ncores = 4, counts = NULL, ids = NULL, verbose = TRUE){
    specific.genes <- names(entropies)

    # only differentially expressed genes with low entropy are considered
    if(any(rownames(differential.genes) %in% specific.genes)){
        differential.genes <- differential.genes[which(rownames(differential.genes) %in% specific.genes),]
    }
    if(verbose) cat(nrow(differential.genes), " genes in intersect of differential and specific genes\n", sep = "")
    # genes in data frame are already in order

    # create two or three gene lists
    # in all lists genes are in same order, but their assigned scores represent
    # their position in the ranking of that score
    # total ranking is the sum of all rankings
    dsg_ranks <- 1:nrow(differential.genes)
    names(dsg_ranks) <- rownames(differential.genes)
    
    specific.ranks <- 1:length(dsg_ranks)
    names(specific.ranks) <- names(sort(entropies[names(dsg_ranks)]))
    specific.ranks <- specific.ranks[names(dsg_ranks)]
    
    # if lasso data is available include lls in the ranking
    # else calculate lasso or do not rank by lls
    if(is.null(lasso.data)){
	    if(!build.lasso){
        	ranks <- deg.weight * dsg_ranks + entropy.weight * specific.ranks
        	names(ranks) <- names(dsg_ranks)
        	ranks <- sort(ranks)
	    }else{
		    if(!is.null(counts) && !is.null(ids)){
			    if(verbose) cat("Calculating lasso models for ", length(dsg_ranks), " genes...\n", sep = "")
		    	lasso.data <- build_lassos(counts[, names(dsg_ranks)], ids, "", NULL, ncores, gamma)
		    }
	    }
    }
    if(!is.null(lasso.data)){
        # extract lls data for genes in the count matrix
        lls <- lasso.data$lls
        lls <- lls[names(dsg_ranks)]
        lls_ranks <- 1:length(lls)
        names(lls_ranks) <- names(lls)[order(lls)]
        lls_ranks <- lls_ranks[names(dsg_ranks)]
        
        ranks <- lls.weight * lls_ranks + deg.weight * dsg_ranks + entropy.weight * specific.ranks
        names(ranks) <- names(dsg_ranks)
        ranks <- sort(ranks)
    }

    if(!is.null(lasso.data)){
        # keep only genes with fitting models
	      if(verbose) cat("Filtering for lls fit...\n")
        differential.genes <- differential.genes[names(lls)[which(lls < summary(lls)[4])],]
	      if(verbose) cat(nrow(differential.genes), " genes passed lls threshold\n", sep = "")
    }
    
    ranks <- ranks[which(names(ranks) %in% rownames(differential.genes))]
    differential.genes <- differential.genes[names(ranks),]
    
    if(!is.null(lasso.data)){
        differential.genes <- cbind(differential.genes, lls = lls[names(ranks)], entropy = entropies[names(ranks)])
    }else{
        differential.genes <- cbind(differential.genes, entropy = entropies[names(ranks)])
    }

    return(differential.genes)
}
