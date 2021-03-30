#' this function is a wrapper for EBImage function readImage
#' read ST image and scale to 1000 x 1000
#'
#' @param path filepath to ST slide image

read_image <- function(path, hw = 1000){
  # check input
	if(!is.character(path)){
		warning("'path' must be a character")
		return(NULL)
	}
	if(!file.exists(path)){
		warning("file could not be found")
		return(NULL)
	}
  # read
	img <- try(EBImage::readImage(path))
	if(class(img) == "try-error"){
		warning("EBImage could not read specified file")
		return(NULL)
	}

	# resize to 1000 x 1000, 
	# enough for good plots and easier on memory
	return(EBImage::resize(img, h = hw, w = hw))
}
