#' fill the barcode data frame with artificial entries for missing coordinates
#'
#' @param ids barcode data frame assigning spatial coordinates to all spots; column names must be 'X' and 'Y'
#' @return full ids data frame

fill_ids <- function(ids) {
    rnames <- c()
    add_ids <- c()
    for(i in 2:34){
        for(j in 2:32){
            idx <- intersect(which(ids[,"X"] == i), which(ids[,"Y"] == j))
            if(length(idx) == 0){
                if(colnames(ids)[1] == "X")
                    add_ids <- rbind(add_ids, c(i, j))
                else 
                    add_ids <- rbind(add_ids, c(j, i))
                rnames <- c(rnames, paste("fill",as.character((i-1)*32+j), sep = "_"))
            }
        }
    }
    if(length(add_ids) > 0){
    	rownames(add_ids) <- rnames
    	ids <- as.data.frame(rbind(as.matrix(ids), as.matrix(add_ids)))
    }
    return(ids)
}
