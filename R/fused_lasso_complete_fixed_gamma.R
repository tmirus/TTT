#' Building sparse fused lasso models for given genes from their count data
#' @details This function will build a sparse fused lasso for each input gene from it's expression vector and the associated spatial coordinates and save the resulting models to an .Rdata object
#' @param counts_matrix a matrix with gene expression levels across genes (columns) and associated spots (rows)
#' @param ids_table table  associating each rowname of the counts_matrix with an x and y coordinate
#' @param name name of the object to save to
#' @param output_folder folder to save object in
#' @param gamma ratio of penalization between sparsity and continuity (default 1, high = sparse)
#' @return list containing three elements\cr
#' 1) counts - matrix containing the fits (Y) for each gene's model
#' 2) BICs - vector containing the BIC for each model
#' 3) lls - vector containing the likelihood for each model

fused_lasso_complete_fixed_gamma <- function(counts_matrix, ids_table, name, output_folder, gamma=1){
    gene_list <- colnames(counts_matrix)

    # store useful information from the genlasso model(s)
    tmp_fits_set <- c()
    #fits_set <- c()
    BICs_set <- c()
    lls_set <- c()

    # not all genes can be used
    to.remove <- c()

    # serial calculation of models for all genes in counts_matrix
    for (k in seq_len(length(gene_list))) {
      i <- gene_list[k]
      # create a 2D matrix (x and y information) for each gene for use with genlasso
      # y in rows, x in columns (by our definition)
      tmp_matrix <- matrix(data = 0,nrow = max(ids_table[,1]),ncol = max(ids_table[,2]))
      for (j in 1:nrow(counts_matrix)){
          tmp_matrix[ids_table[rownames(counts_matrix)[j],1],ids_table[rownames(counts_matrix)[j],2]] <- counts_matrix[j,i]
      }

      # full solution path calculation
      tmp_lasso <- fusedlasso2d(tmp_matrix,gamma = gamma)
      
      # calculate BICs for all possible models
      BIC_list <- calc_BICs(tmp_lasso)

      # if there is no model with minimal BIC
      # the model building failed; skip those genes
      if(length(which.min(BIC_list[,4])) > 0){
      	# retrieve lambda corresponding to minimal BIC
        lambda_BIC <- BIC_list[which.min(BIC_list[,4]),2]
        
        # retrieve the model corresponding to the lambda
        # coef_error catches truncation errors and 
        # uses a different but similar lambda
        fits <- coef_error(tmp_lasso, lambda=lambda_BIC)
        tmp_fit <- fits$beta

        # contains the count models (Y)
        tmp_fits_set <- cbind(tmp_fits_set,tmp_fit)
        
        # contains the fit objects
        #fits_set <- append(fits_set,list(fits))
        
        # contains the BICs
        BICs_set <- append(BICs_set,min(BIC_list[,4]))
        
        # contains the model fit ratio (part of BIC)
        lls_set <- append(lls_set,min(BIC_list[,5]))
      }else{
  	to.remove <- c(to.remove, which(gene_list == i))
      }
    }
    # if any genes failed, remove their names
    if(length(to.remove)>0){
      # if all genes are removed, return empty list
      if(length(to.remove) < length(gene_list)){
        colnames(tmp_fits_set) <- gene_list[-to.remove]
        names(lls_set) <- gene_list[-to.remove]
      }else{
        return(list())
      }
    }else{
      colnames(tmp_fits_set) <- gene_list
      names(lls_set) <- gene_list
    }
    return(list(counts = tmp_fits_set, #fits = fits_set, 
                BICs = BICs_set, 
                lls = lls_set)
           )
}
