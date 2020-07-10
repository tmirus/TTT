#' this function returns a table of interesting genes based on lasso model fit, specificity and differential expression
#' 
#' @param analysis.data output of analyze_clustering function
#' @param lasso.data output of build_lassos function
#' @return data frame containing gene name, cluster with lowest p-value, lowest p-value and up/downregulation information for each gene
#' @export
filter_genes <- function(analysis.data, lasso.data){
    differential.genes <- analysis.data$differential_genes
    specific.genes <- names(analysis.data$specific)
    lls <- lasso.data$lls

    # only differentially expressed genes with low entropy are considered

    # maybe reorder by combined rank? combined score?
    if(any(rownames(differential.genes) %in% specific.genes)){
        differential.genes <- differential.genes[which(rownames(differential.genes) %in% specific.genes),]
    }

    # plot lls gene ranks vs combined gene ranks

    #plot("test_lls_dist.pdf")
    #hist(lls, breaks = 200, main = "lls all genes")
    #plot(density(lls), col = "red")
    #abline(v = mean(lls), col = "red")
    lls <- lls[rownames(differential.genes)]
    lls <- sort(lls)
    #lines(density(lls), col = "green")
    #abline(v = mean(lls), col = "green")
    
    lls_ranks <- 1:length(lls)
    names(lls_ranks) <- names(lls)
    
    dsg_ranks <- 1:length(lls)
    names(dsg_ranks) <- rownames(differential.genes)
    dsg_ranks <- dsg_ranks[names(lls)]

    specific.ranks <- 1:length(lls)
    names(specific.ranks) <- names(sort(analysis.data$specific[names(lls)]))
    specific.ranks <- specific.ranks[names(lls)]

    ranks <- lls_ranks + dsg_ranks + specific.ranks
    names(ranks) <- names(lls)
    ranks <- sort(ranks)

    #hist(lls, breaks = 100, main = "lls dsg")
    #plot(lls_ranks, dsg_ranks, xlab = "rank(LLS)", ylab = "rank(DSG)", main = "relation between LLS and DSG")
    #dev.off()

    differential.genes <- differential.genes[names(lls)[which(lls < summary(lls)[2])],]
    ranks <- ranks[which(names(ranks) %in% rownames(differential.genes))]
    differential.genes <- differential.genes[names(ranks),]

    return(differential.genes)
}
