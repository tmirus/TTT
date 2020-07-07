
fill_ids <- function(ids) {
    rnames <- c()
    add_ids <- c()
    for(i in 2:34){
        for(j in 2:32){
            idx <- intersect(which(ids[,"X"] == i), which(ids[,"Y"] == j))
            if(length(idx) == 0){
                add_ids <- rbind(add_ids, c(i, j))
                rnames <- c(rnames, paste("fill",as.character((i-1)*34+j), sep = "_"))
            }
        }
    }
    rownames(add_ids) <- rnames
    ids <- as.data.frame(rbind(as.matrix(ids), as.matrix(add_ids)))
    return(ids)
}