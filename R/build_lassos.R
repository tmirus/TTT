#' Calculate 2D fused-lasso solutionpath for a countmatrix and save as object
#' @details This function builds a solution path for each gene (column) in a input matrix, selects the models according to the BIC and 
#' builds a count matrix
#' @param counts input matrix of counts with genes in columns and spot identifier in rows
#' @param ids coordinate table which converts the spot identifier into an x and y coordinate on a grid
#' @param name Name of the experiment
#' @param output_folder folder in which to save the results
#' @param ncores number of parallel processes to start (default 8)
#' @param gamma ratio of penalization between sparsity and continuity (default 1, high = sparse)
#' @return list containing two elements\cr
#' 1) counts - count matrix in the same format as input count matrix
#' 2) ids - id data frame that assigns spot identifiers to locations on the grid
#' @export

build_lassos <- function(counts,ids,name,output_folder = NULL,ncores=4,gamma=1){
    suppressMessages(library(parallel, quietly = TRUE))
    # select genes randomly for each thread to balance the load
    genes <- colnames(counts)
    gene_list <- vector(mode = "list", length = ncores)


    non.zero <- sum(rowSums(counts) > 0)
    zero.fraction <- sum(counts == 0)/length(as.vector(counts))

    if(ncores > 1){
        genes.per.thread <- as.integer(ncol(counts) / ncores) 
        for(i in 1:(ncores-1)){
            indices <- sample(length(genes), genes.per.thread, replace = FALSE)
            gene_list[[i]] <- genes[indices]
            genes <- genes[-indices]
        }
    }
    gene_list[[ncores]] <- genes

    # building the lasso model with ncores processes in parallel
    lasso.data <- mclapply(seq_len(ncores),
                           function(x) 
                               fused_lasso_complete_fixed_gamma(
                                   counts = counts[, gene_list[[x]], drop=F],
                                   ids_table = ids,
                                   name = paste(name,x,sep='_'),
                                   output_folder = output_folder,
                                   gamma=gamma
                               ),
                           mc.cores = ncores
                           )
    
    # reconstruct the counts matrix and store additional information 
    # returned by fused_lasso_complete_fixed_gamma
    counts <- c()
    fits <- c()
    BICs <- c()
    lls <- c()
    for(i in 1:length(lasso.data)){
        if(length(lasso.data) > 0){
            counts <- cbind(counts, lasso.data[[i]]$counts)
            #fits <- append(fits, lasso.data[[i]]$fits)
            BICs <- append(BICs, lasso.data[[i]]$BICs)
            lls <- append(lls, lasso.data[[i]]$lls)
        }
    }
    
    if (!is.null(output_folder) && is.character(output_folder)) {
        # create output folder if it does not exist
        if(!dir.exists(paste(output_folder, name, sep = '/'))){
            dir.create(paste(output_folder, name, sep = '/'), recursive = TRUE)
        }
        
        # save the counts matrix and additional information to output folder
        saveRDS(as.matrix(counts), 
                file = paste(output_folder, name, "counts_lasso_BIC.RDS", sep = '/'))
        saveRDS(fits, file = paste(output_folder, name, "fits_lasso_BIC.RDS", sep ='/'))
        saveRDS(BICs, file = paste(output_folder, name, "BICs_lasso_BIC.RDS", sep = '/'))
        saveRDS(lls, file = paste(output_folder, name, "lls_lasso_BIC.RDS", sep = '/'))
    }

    # keep the original ids and rownames
    # create a lasso ids data frame based on mapping of index to x and y
    ids_lasso <- c()
    for(i in 0:(max(ids[,2])-1)){
        for(j in 1:max(ids[,1])){
            ids_lasso <- rbind(ids_lasso, c(j, i))
        }
    }

    rownames(ids_lasso) <- 1:nrow(ids_lasso)
    ids_lasso[,2] <- ids_lasso[,2] + 1
    # match the temporary data frame to original ids
    lasso.data <- prepare_lasso_data(ids, ids_lasso, as.matrix(counts))
    lasso.data[["lls"]] <- lls
    
    rm(counts)
    
    if (!is.null(output_folder) && is.character(output_folder)) {
        saveRDS(as.matrix(lasso.data$counts), 
                file = paste(output_folder, name, "counts_lasso_BIC.RDS", sep = '/'))
        saveRDS(lasso.data$ids,
                file = paste(output_folder, name, "ids.RDS", sep = '/'))
    }
    
    non.zero.after <- sum(rowSums(lasso.data$counts) > 0)
    if(non.zero.after < 0.75 * non.zero){
        warning("More than 25% of non-zero spots have been set to 0 by lasso. Consider using a smaller gamma for less sparsity.")
    }
    
    zero.fraction.after <- sum(lasso.data$counts == 0)/length(as.vector(lasso.data$counts))
    cat("Fraction of zeros in raw data: ", zero.fraction, "\n", sep = "")
    cat("Fraction of zeros in lasso data: ", zero.fraction.after, "\n", sep = "")
    
    zero.genes <- which(colSums(lasso.data$counts) == 0)
    if(length(zero.genes) > 0){
        cat(length(zero.genes), " genes have been set to 0. Removing...\n", sep = "")
        lasso.data$counts <- lasso.data$counts[, -zero.genes]
        lasso.data$lls <- lasso.data$lls[colnames(lasso.data$counts)]
    }
    zero.spots <- which(rowSums(lasso.data$counts) == 0)
    if(length(zero.spots) > 0){
        lasso.data$counts <- lasso.data$counts[-zero.spots, ]
    }

    # return counts and ids
    return(lasso.data)
}
