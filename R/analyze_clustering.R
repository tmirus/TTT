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
    cluster.libs <- cluster.libs[, which(colSums(cluster.libs) > 0)]
    cluster.libs <- apply(cluster.libs, 2, function(x){x / sum(x)}) 

    specific <- sapply(colnames(cluster.libs), function(x){
        if(! sum(cluster.libs[,x]) > 0) return(NA)
        if(any(cluster.libs[,x] == 0)){
            entropy <- -sum(cluster.libs[-which(cluster.libs[,x] == 0), x] * log2(cluster.libs[-which(cluster.libs[,x] == 0), x]))
        }else{
            entropy <- -sum(cluster.libs[, x] * log2(cluster.libs[, x]))
        }
        return(entropy)
    })

    names(specific) <- colnames(cluster.libs)
    print(str(colnames(cluster.libs)))
    if(any(is.na(specific))){
    	specific <- specific[-which(is.na(specific))]
    }
    specific <- specific[which(specific <= 1)]
    #pdf("/home/tmirus/TTT/development/entropy_test.pdf")
    #plot(density(specific))
    #hist(specific, breaks = 150)
    #dev.off()
    #print(summary(specific))

    # implement testing with multtest for differentially expressed genes for each cluster
    test.results <- list()
    for(cl1 in rownames(cluster.libs)){
        temp.counts <- t(counts[c(which(clustering == cl1),which(clustering != cl1)),])
        labs <- c(rep(0,length(which(clustering == cl1))), rep(1,length(which(clustering != cl1))))
        t.scores <- mt.teststat(temp.counts, labs, test = "t")
        p.vals <- 2 * pt(abs(t.scores), length(labs) - 2, lower.tail = FALSE)
        p.vals <- p.adjust(p.vals, method = "BH")
        names(p.vals) <- rownames(temp.counts)
        test.results[[cl1]] <- p.vals
    }
    saveRDS(test.results, "test_results.RDS")

    gene.lists <- list()
    for(cl1 in rownames(cluster.libs)){
        gene.scores <- test.results[[cl1]]
        names(gene.scores) <- names(test.results[[cl1]])

	print(length(gene.scores))
	print(sum(is.na(gene.scores)))
        if(any(is.na(gene.scores))){
            gene.scores[is.na(gene.scores)] <- 1
        }
	print(sum(is.nan(gene.scores)))
        if(any(is.nan(gene.scores))){
            gene.scores[is.nan(gene.scores)] <- 1
        }
	print(sum(is.na(gene.scores)))
        gene.lists[[cl1]] <- sort(gene.scores)
    }
	saveRDS(gene.lists, "gene_lists.RDS")
    all.genes <- as.vector(sapply(gene.lists, function(x){names(x)}))
    print(str(all.genes))
    all.genes <- all.genes[!is.na(all.genes)]
    all.genes <- unique(all.genes)
    print(str(all.genes))
    saveRDS(all.genes, "all_genes.RDS")

    clusters <- sort(unique(clustering))
    print(clusters)
    gene.table <- c()
    for(g in all.genes){
        cluster.values <- c()
        for(cl in as.character(clusters)){
		#print(str(names(gene.lists[[cl]])))
            if(g %in% names(gene.lists[[cl]])){
                cluster.values <- c(cluster.values, as.numeric(gene.lists[[cl]][g]))
            }else{
                cluster.values <- c(cluster.values, 1)
            }
        }
        gene.table <- rbind(gene.table, cluster.values)
    }
    saveRDS(gene.table, "gene_table.RDS")
    print("gene.table")
    print(dim(gene.table))
    gene.table <- as.data.frame(gene.table)
    colnames(gene.table) <- clusters
    rownames(gene.table) <- all.genes
    for(cl in colnames(gene.table)){
	    gene.table[[cl]] <- as.numeric(as.character(gene.table[[cl]]))
    }
	print("set rownames")
    to.remove <- c()
    for(i in 1:nrow(gene.table)){
        if(min(gene.table[i,]) > sig.level | is.na(min(gene.table[i,]))){
            to.remove <- c(to.remove, i)
        }
    }
    print("removing rows")
    if(length(to.remove) > 0)
        gene.table <- gene.table[-to.remove,]
    print(dim(gene.table))
    print("done")
    gene.table <- gene.table[order(apply(as.matrix(gene.table), 1, min)),]
    #gene.table <- gene.table[-which(is.na(apply(as.matrix(gene.table),1,min))),]
    print(str(gene.table))

    print("working until here")
    saveRDS(gene.table, "gene_table2.RDS")
    gene.info <- c()
    for(i in 1:nrow(gene.table)){
        g <- rownames(gene.table)[i]
        p <- min(gene.table[i,])
        cl <- clusters[which.min(gene.table[i,])]
	cat(g, p, cl, sep = "\t")
	if(length(cl) > 0){
        	reg <- sign(mean(counts[which(clustering == cl),g]) - mean(counts[which(clustering != cl), g]))
		cat("\t",reg, "\n")
        	gene.info <- rbind(gene.info,
                	           c(g, cl, p, reg))
	}
    }
    gene.info <- as.data.frame(gene.info)
    print("gene.info")
    print(str(gene.info))
    colnames(gene.info) <- c("gene", "cluster", "pVal", "regulation")
    gene.info$pVal <- as.numeric(as.character(gene.info$pVal))
    gene.info$regulation <- as.numeric(as.character(gene.info$regulation))
    print(str(gene.info$gene))
    rownames(gene.info) <- as.character(gene.info$gene)
    print("done")

    return(list(specific_genes = names(specific), differential_genes = gene.info, dsg = gene.info[which(rownames(gene.info) %in% names(specific)),]
                ))
}
