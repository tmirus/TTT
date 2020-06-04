#' @export
analyze_clustering <- function(counts, ids, clustering.info, specificity.threshold = 0.9, sig.level = 0.0001){
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

    # implement testing with multtest for differentially expressed genes for each cluster
    test.results <- list()
    for(cl1 in rownames(cluster.libs)){
        test.results[[cl1]] <- list()
        for(cl2 in rownames(cluster.libs)){
            if(cl1 == cl2) next
            temp.counts <- t(counts[c(which(clustering == cl1),which(clustering == cl2)),])
            labs <- c(rep(0,length(which(clustering == cl1))), rep(1,length(which(clustering == cl2))))
            
            t.scores <- mt.teststat(temp.counts, labs, test = "t")
            p.vals <- 2 * pt(abs(t.scores), length(labs) - 2, lower.tail = FALSE)
            p.vals <- p.adjust(p.vals, method = "BH")
            names(p.vals) <- colnames(cluster.libs)
            test.results[[cl1]][[cl2]] <- p.vals
        }
    }

    gene.lists <- list()
    for(cl1 in rownames(cluster.libs)){
        genes <- c()
        for(cl2 in rownames(cluster.libs)){
            genes <- c(genes, names(test.results[[cl1]][[cl2]])[which(test.results[[cl1]][[cl2]] < sig.level)])
        }
        genes <- unique(genes)
        gene.mat <- matrix(1, nrow = nrow(cluster.libs)-1, ncol = length(genes))
        rownames(gene.mat) <- rownames(cluster.libs)[-which(rownames(cluster.libs) == cl1)]
        colnames(gene.mat) <- genes
        for(cl2 in rownames(gene.mat)){
            gene.mat[cl2,] <- test.results[[cl1]][[cl2]][genes]
        }
        if(any(is.na(gene.mat))){
            gene.mat[is.na(gene.mat)] <- 1
        }
        if(any(is.nan(gene.mat))){
            gene.mat[is.nan(gene.mat)] <- 1
        }
        gene.scores <- colMeans(gene.mat)
        gene.lists[[cl1]] <- sort(gene.scores)
    }

    return(list(specific, gene.lists))
}