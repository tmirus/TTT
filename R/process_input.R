#' process input data from ST pipeline to make sure they conform 
#' to the format the package expects
#' 
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param ids data frame or matrix assigning spatial coordinates to the spots,
#' the spot 0,0 is in the right lower corner of the 
#' corresponding tissue slide and the number of spots along the X-axis is 
#' higher than the number of spots on the Y-axis.
#' optional if spatial information is contained in rownames of counts
#' @param separate_by only necessary if spatial information is given as
#' in the output of ST Pipeline (Navarro et al.), e.g. 32x2; default x
#' @param force_counts logical, in case there are less genes than spots this needs to be TRUE; 
#' default FALSE
#' @param force_indices logical, in case the range of x and y indices is not the same as in the
#' ST barcode files
#' @param dup.sep character, optional; if duplicate gene names are distinguished by e.g. _1, _2
#' at the end, pass dup.sep="_" to try to combine these columns
#' @param n.gene.cutoff integer > 0, number of genes that need to be expressed in a spots in order 
#' for the spot not to be removed from the data
#' @return list with two entries:\cr
#' 1) counts - count matrix
#' 2) ids - barcode data frame assigning spatial positions to spots
#' @export

process_input <- function(counts, ids = NULL, separate_by = "x", force_counts = FALSE, force_indices = FALSE, dup.sep = NULL, n.gene.cutoff = 100){
  if(is.character(counts)){
    if(file.exists(counts)){
      counts <- as.matrix(read.table(counts, check.names = FALSE))
      if(!is.numeric(counts)){
        stop("An error occurred while reading the counts matrix. Try reading it manually and pass it as parameter.")
      }
    }
  }
  # check input parameters
  if(!is.matrix(counts)){
    stop("counts needs to be supplied as matrix")
  }
  if(!is.numeric(counts) || any(counts < 0)){
    stop("counts must be a non-negative numeric matrix")
  }
  if(nrow(counts) >= ncol(counts) && !force_counts){
    stop("There seem to be more spots than genes. Please check matrix
         orientation or set force_counts = TRUE to use the data as it is.")
  }
  if(n.gene.cutoff <= 0 || round(n.gene.cutoff) != n.gene.cutoff){
    stop("Invalid value for 'n.gene.cutoff'. Must be integer > 0.")
  }
  if(!is.null(ids)){
    if(! (is.matrix(ids) || is.data.frame(ids)) ){
      stop("ids must be a matrix or a data frame")
    }
    if(ncol(ids) != 2){
      stop("ids must contain two columns (X and Y)")
    }
    # fill ids to range from 2 to 32/34 in Y/X-direction
    ids <- fill_ids(ids)
    
    if(! all(range(ids) == c(2,34) && !force_indices) ){
      stop("The range of indices in ids does not match typical ST barcode files (index range 2 - 34). 
           Set force_indices = TRUE to proceed anyways. This might lead to unexpected behaviour.")
    }
    if(nrow(ids) < nrow(counts)){
      stop("ids does not have enough rows to contain spatial information for all spots in counts.")
    }
    if(!all(rownames(counts) %in% rownames(ids))){
      stop("One or more spots are missing spatial information.")
    }
  }else{
    if(! all(grep(pattern=separate_by, x = rownames(counts)) == 1:nrow(counts))){
      stop("Splitting character split_by does not seem to be present in all spot names.")
    }
  }
  
  # prepare either input ids or convert rownames to ids data frame
  if(!is.null(ids)){
    if(!all(colnames(ids) == c("X", "Y"))){ 
      # reorder the ids data frame to have x-coordinates in second and
      # y-coordinate in first column
      input.ids <- ids
      r1 <- range(ids[,1])
      r2 <- range(ids[,2])
      if( (r1[2] - r1[1]) > (r2[2] - r2[1]) ){
        x_ind <- 1
        y_ind <- 2
      }else if( (r1[2] - r1[1]) < (r2[2] - r2[1]) ){
        x_ind <- 2
        y_ind <- 1
      }else{
        stop("Dimension of a ST slide is 34x32. X and Y in the id matrix have the same dimension.")
      }
      ids <- data.frame(Y = input.ids[, y_ind], X = input.ids[, x_ind])
      rownames(ids) <- rownames(input.ids)
    }
  }else{
    # extract ids from rownames of counts and create a data frame
    ids.parts <- strsplit(rownames(counts), split = separate_by)
    coords1 <- as.numeric(sapply(ids.parts, function(x){x[[1]]}))
    coords2 <- as.numeric(sapply(ids.parts, function(x){x[[2]]}))
    
    if(any(is.na(c(coords1,coords2)))){
      stop("rownames of counts could not be converted to coordinates. Make sure the rownames
           conform to the pattern {coord1}{split_by}{coord2} or provide an ids matrix
           mapping barcodes to spatial positions")
    }
    if(any(c(coords1, coords2) < 0)){
      stop("Invalid ids extracted from rownames. Coordinates must be >= 0.")
    }
    
    r1 <- range(coords1)
    r2 <- range(coords2)
    if( (r1[2] - r1[1]) > (r2[2] - r2[1]) ){
      ids <- data.frame(Y = coords2, X = coords1)
    }else if( (r1[2] - r1[1]) < (r2[2] - r2[1]) ){
      ids <- data.frame(Y = coords1, X = coords2)
    }else{
      stop("Dimension of a ST slide is 34x32. X and Y in the id matrix have the same dimension.")
    }
    rownames(ids) <- rownames(counts)
  }
  
  # adjust coordinates for compatibility 
  # with the package if possible
  ids[,"X"] <- ids[,"X"] + (2 - min(ids[,"X"]))
  ids[,"Y"] <- ids[,"Y"] + (2 - min(ids[,"X"]))
  
  if(!all(range(ids[,1]) == c(2, 32)) || !all(range(ids[,2]) == c(2,34)) ){
    stop("Could not adjust data to conform to package standards. Please check your input.")
  }

  # deal with potential duplicates
  if(!is.null(dup.sep)){
    gene.names <- strsplit(colnames(counts), split = dup.sep)
    suppressWarnings(gene.names <- sapply(gene.names, function(x){
    if(!is.na(as.numeric(x[length(x)])) && !(nchar(x[length(x)]) > 1)){
        if(length(x) == 2){
        return(x[1])
      }else{
        return(paste(x[-length(x)], collapse = "_"))
      }
    }
    return(paste(x, collapse = "_"))
    }))
    colnames(counts) <- gene.names
  }
  cat("Duplicates found:", sum(duplicated(colnames(counts))), "\nRemoving...\n")
  
  # find all duplicate names, iterate over them, sum them up and remove superfluous columns
  if(sum(duplicated(colnames(counts))) > 0){
    dup.genes <- unique(colnames(counts)[duplicated(colnames(counts))])
    for(g in dup.genes){
      gene.idx <- which(colnames(counts) == g)
      counts[,gene.idx[1]] <- rowSums(counts[, gene.idx])
      counts <- counts[,-gene.idx[-1]]
    }
  }
  
  # remove spots with exceptionally low number of expressed genes
  n.genes <- apply(counts, 1, function(x) sum(x>0))
  if(any(n.genes) <= 100){
    counts <- counts[-which(n.genes < 100),]
  }
  
  return(list(counts = counts, ids = ids))
}
