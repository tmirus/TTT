#' process the lasso counts and ids for conformity with pipeline output
#' 
#' @description after model building the resulting counts matrix and the 
#' ids data frame do not contain rownames (spot names / spatial information).
#' This function matches spatial positions of lasso and raw ids and 
#' integrates the information into the lasso model.
#' @param ids_raw data frame containing spatial information as built by
#' load_ids from and ST barcode file
#' @param ids_lasso data frame containing ids of lasso spots. For this the 
#' order in which spots are stored in the lasso models needs to be known:
#' idx = i*x_max + j; 0 <= i <= x_max-1, 1 <= j <= y_max
#' @param counts_lasso the model that is built in build_lassos, 
#' containing 34*32 spots without rownames
#' @return list containing two elements\cr
#' 1) counts - the count matrix with rownames and reduced to original spots
#' 2) ids - the lasso ids

prepare_lasso_data <- function(ids_raw, ids_lasso, counts_lasso){
    # store the indices corresponding to matching spots
    idx_lasso <- c()
    idx_raw <- c()

    # iterate over spots in lasso ids and search matching rows in raw ids
    for(i in 1:nrow(ids_lasso)){
        spot.intersect <- intersect(
            which(ids_raw[,"Y"] == ids_lasso[i,1]), 
            which(ids_raw[,"X"] == ids_lasso[i,2])
        )
        if(length(spot.intersect) == 1){
            idx_lasso <- c(idx_lasso, i)
            idx_raw <- c(idx_raw, spot.intersect)
        }
    }
    # reduce the lasso ids and counts to spots that are contained in raw ids
    ids_lasso <- ids_lasso[idx_lasso,]
    counts_lasso <- counts_lasso[idx_lasso,]
    
    # name the spots accordingly
    rownames(ids_lasso) <- rownames(ids_raw)[idx_raw]
    rownames(counts_lasso) <- rownames(ids_raw)[idx_raw]
    colnames(ids_lasso) <- c("Y", "X")
    
    return(list(counts=counts_lasso, ids = ids_lasso))
}
