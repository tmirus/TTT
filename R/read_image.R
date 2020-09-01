#' this function is a wrapper for EBImage function readImage
#'
#' @param path filepath to ST slide image

read_image <- function(path){
	if(!is.character(path)){
		warning("'path' must be a character")
		return(NULL)
	}
	if(!file.exists(path)){
		warning("file could not be found")
		return(NULL)
	}
	img <- try(EBImage::readImage(path))
	if(class(img) == "try-error"){
		warning("EBImage could not read specified file")
		return(NULL)
	}

	return(EBImage::resize(img, h = 1000, w = 1000))
}
