#' remove background spots (not covered by tissue) from ST data based on tissue image.
#' Manual checking of the result and adjustment of the threshold is needed.
#' 
#' @param img image of the ST slide, cropped to the spot area. image read by EBimage::readImage()
#' @param nx number of spots shown in the image on the horizontal axis
#' @param ny number of spots shown in the image on the vertical axis
#' @param ids data frame assigning barcodes / spot names to spatial coordinates
#' @param counts non-negative numeric matrix containing gene counts, 
#' rows correspond to spots, columns correspond to genes
#' @param threshold relative brightness. spots above this threshold are discarded
#' @details nx and ny depend on the cropping of the image. the area of measurements in an ST experiment is
#' @return list containing two entries\cr
#' 1) spots.to.keep - vector of barcodes / spot names of the spots that should not be removed
#' 2) spots.keep.clustering - vector of barcodes / spot names of the spots that should not be removed based on clustering in addition to image analysis
#' 3) image - image showing the areas that are retained.
#' typically framed by spots for which there exist no measurements. Depending on whether these are contained 
#' in the image or not, different values need to be chosen for nx, ny. If they are included, the 
#' default values of nx=35 and ny=33 are sensible. If not, it should be nx=33 and ny=31.\cr
#' The image is assumed to be oriented such that the rectangular 4x4 array of spots in one of the corners
#' is in the lower right corner and on the x-axis (left to right) there are more spots than on the y-axis.
#' @export

remove_background <- function(img, nx = 35, ny=33, ids, counts, threshold = 0.7){
  # reduce to one pixel per spot
  img <- EBImage::channel(img, "gray")
  img.reduced <- EBImage::resize(img, w=nx,h=ny)
  names <- rownames(counts)
  
  # blur to reduce the risk of removing spots within tissue
  for(x in 1:nx){
    for(y in 1:ny){
      x.neighbours <- seq(x-1,x+1)
      y.neighbours <- seq(y-1,y+1)
      if(any(x.neighbours < 1 | x.neighbours > nx)) 
        x.neighbours <- x.neighbours[-which(x.neighbours < 1 | x.neighbours > nx)]
      if(any(y.neighbours < 1 | y.neighbours > ny)) 
        y.neighbours <- y.neighbours[-which(y.neighbours < 1 | y.neighbours > ny)]
      imageData(img.reduced)[x,y] <- mean(imageData(img.reduced)[x.neighbours, y.neighbours])
    }
  }
  
  # select the spots to remove
  imageData(img.reduced)[imageData(img.reduced) > threshold*max(imageData(img.reduced))] <- 1
  to.remove <- t(matrix(imageData(img.reduced) == 1, nrow=nx))
  
  # store in more compatible format
  spots.to.keep <- c()
  for(x in 1:nx){
    for(y in 1:ny){
      if(!to.remove[y,x]) spots.to.keep <- rbind(spots.to.keep, c(35-x,33-y))
    }
  }
  spots.to.keep <- as.data.frame(spots.to.keep)
  colnames(spots.to.keep) <- c("X", "Y")
  
  img <- EBImage::resize(img, 1000, 1000)
  # return the original image with blanks where spots were removed
  dx <- dim(img)[1] / nx
  dy <- dim(img)[2] / ny
  for(x in 1:nx){
    for(y in 1:ny){
      if(to.remove[y,x]){
        imageData(img)[(as.integer((x-1)*dx+1)):as.integer(x*dx), (as.integer((y-1)*dy+1)):as.integer(y*dy)] <- 1
      }
    }
  }
  if(!is.null(ids) && !is.null(names)){
    spots <- sapply(names,
                    function(x){
                      if(length(intersect(which(spots.to.keep[,"X"] == ids[x,2]), which(spots.to.keep[,"Y"] == ids[x,1])))>0){
                        return(TRUE)
                      }else{
                        return(FALSE)
                      }
                    }
    )
    spots <- names(spots)[spots]
    spots.to.keep <- spots
  }
  ids.reduced <- ids[spots,]
  ids.reduced <- clean_ids(ids.reduced)
  spots.to.keep <- rownames(ids.reduced)

  # tSNE embedding
  tsne <- Rtsne(counts, check_duplicates = F)$Y
  # create comparatively large number of clusters
  clustering <- kmeans(tsne, 15)$cluster
  names(clustering) <- rownames(counts)

  # keep all clusters with a certain amount of spots in spots.to.keep
  clusters.to.keep <- c()
  for(cl in unique(clustering)){
    if(sum(names(clustering[clustering == cl]) %in% spots.to.keep) > 0.75 * length(which(clustering == cl))){
      clusters.to.keep <- c(clusters.to.keep, names(clustering[clustering == cl]))
    }
  }

  ids.reduced <- ids[clusters.to.keep,]
  ids.reduced <- clean_ids(ids.reduced)
  clusters.to.keep <- rownames(ids.reduced)

  return(list(
	spots.to.keep = spots.to.keep, 
	spots.keep.clustering = clusters.to.keep, 
	image = img
	))
}
