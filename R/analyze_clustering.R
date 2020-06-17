#' @export
analyze_clustering <- function(counts, ids, clustering.info, specificity.threshold = 0.9, sig.level = 0.05){
    counts <- counts[rownames(clustering.info$scores),]
    clustering <- clustering.info$clustering

    # find cluster-specific genes by calculating total RNA per cluster and gene
    cluster.libs <- matrix(0, ncol=ncol(counts), nrow = length(unique(clustering)))
    colnames(cluster.libs) <- colnames(counts)
    rownames(cluster.libs) <- unique(clustering)

    for(cl in rownames(cluster.libs)){
        cluster.libs[cl,] <- colSums(counts[which(clustering == as.numeric(cl)),,drop=F])
    }

    specific <- sapply(colnames(cluster.libs), function(x){
        if(! sum(cluster.libs[,x]) > 0) return(-1)
        for(i in 1:nrow(cluster.libs)){
            if(cluster.libs[i,x] / sum(cluster.libs[,x]) > specificity.threshold){
                return(as.numeric(rownames(cluster.libs)[i]))
            }
        }
        return(-1)
    })

    names(specific) <- colnames(cluster.libs)
    specific <- specific[-which(specific < 0)]

	#counts <- log2(counts + 1)

    # implement testing with multtest for differentially expressed genes for each cluster
    test.results <- list()
    for(cl1 in rownames(cluster.libs)){
        temp.counts <- t(counts[c(which(clustering == cl1),which(clustering != cl1)),])
        labs <- c(rep(0,length(which(clustering == cl1))), rep(1,length(which(clustering != cl1))))
        t.scores <- mt.teststat(temp.counts, labs, test = "t")
        p.vals <- 2 * pt(abs(t.scores), length(labs) - 2, lower.tail = FALSE)
        p.vals <- p.adjust(p.vals, method = "BH")
        names(p.vals) <- colnames(cluster.libs)
        test.results[[cl1]] <- p.vals
    }

    gene.lists <- list()
    for(cl1 in rownames(cluster.libs)){
        gene.scores <- test.results[[cl1]][which(test.results[[cl1]] < as.numeric(sig.level))]
        names(gene.scores) <- names(test.results[[cl1]])[which(test.results[[cl1]] < as.numeric(sig.level))]

        if(any(is.na(gene.scores))){
            gene.scores[is.na(gene.scores)] <- 1
        }
        if(any(is.nan(gene.scores))){
            gene.scores[is.nan(gene.scores)] <- 1
        }
        gene.lists[[cl1]] <- sort(gene.scores)
    }

    all.genes <- unique(unlist(sapply(gene.lists, function(x){names(x)}), use.names = F))
    clusters <- sort(unique(clustering))
    gene.table <- foreach(g = all.genes, .combine = 'rbind') %do% {
        cluster.values <- c()
        for(cl in clusters){
            if(g %in% names(gene.lists[[cl]])){
                cluster.values <- c(cluster.values, gene.lists[[cl]][g])
            }else{
                cluster.values <- c(cluster.values, 1)
            }
        }
        return(cluster.values)
    }
    gene.table <- as.data.frame(gene.table)
    colnames(gene.table) <- clusters
    rownames(gene.table) <- all.genes

    to.remove <- c()
    for(i in 1:nrow(gene.table)){
        if(min(gene.table[i,]) > sig.level){
            to.remove <- c(to.remove, i)
        }
    }
    if(length(to.remove) > 0)
        gene.table <- gene.table[-to.remove,]


    return(list(specific, gene.lists, gene.table))
}
