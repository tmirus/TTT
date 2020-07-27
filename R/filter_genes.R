#' this function returns a table of interesting genes based on lasso model fit, specificity and differential expression
#' 
#' @param analysis.data output of analyze_clustering function
#' @param lasso.data output of build_lassos function
#' @return data frame containing gene name, cluster with lowest p-value, lowest p-value and up/downregulation information for each gene
#' @export
filter_genes <- function(deg.table, entropies, lasso.data = NULL){
    differential.genes <- deg.table
    specific.genes <- names(entropies)

    # only differentially expressed genes with low entropy are considered
    if(any(rownames(differential.genes) %in% specific.genes)){
        differential.genes <- differential.genes[which(rownames(differential.genes) %in% specific.genes),]
    }

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
    if(!is.null(lasso.data)){
        # extract lls data for genes in the count matrix
        lls <- lasso.data$lls
        lls <- lls[names(dsg_ranks)]
        lls_ranks <- 1:length(lls)
        names(lls_ranks) <- names(lls)[order(lls)]
        lls_ranks <- lls_ranks[names(dsg_ranks)]

        
        ranks <- lls_ranks + dsg_ranks + specific.ranks
        names(ranks) <- names(dsg_ranks)
        ranks <- sort(ranks)
    }else{
        ranks <- dsg_ranks + specific.ranks
        names(ranks) <- names(dsg_ranks)
        ranks <- sort(ranks)
    }


    if(!is.null(lasso.data)){
        # keep only genes with fitting models
        differential.genes <- differential.genes[names(lls)[which(lls < summary(lls)[2])],]
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
