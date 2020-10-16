#' This function helps choosing offset and spacing parameters for good spatial plots. 
#' Basically just a wrapper for spatial_plot that returns the input parameters in a list for easy use later on.
#' 
#' @param ids data frame or matrix assigning spatial coordinates to the spots (should be complete ids data frame, no missing spots)
#' @param img Image of the ST slide, loaded by EBImage::readImage()
#' @param nx numeric > 0, defines spacing between plotted spots on x-axis; 
#' rule of thumb: number of spots shown in the image (horizontal); the large this number, the closer together the spots are 
#' @param ny numeric > 0, defines spacing between plotted spots on y-axis;
#' rule of thumb: number of spots shown in the image (vertical); the larger this number, the closer together the spots are
#' @param ox numeric, defines the offset of spot 0,0 from the lower right corner of the image in x-direction (can be negative)
#' @param oy numeric, defines the offset of spot 0,0 from the lower tight corner of the image in y-direction (can be negative)
#' @param output logical, should an overlay of the adjusted spots be plotted? default TRUE
#' @return list containing the input parameters nx, ny, ox, oy. If the produced plot looks good, this can be used to specify the parameters for plotting in other functions of this package.
#' @export
plot_adjustment <- function(ids, img, nx = 35, ny = 33, ox = 0, oy = 0, output = TRUE){
	if(output){
    plot(
	 spatial_plot(
		rownames(ids), 
		ids, 
		rep(1, nrow(ids)), 
		img, 
		"discrete", 
		list(nx = nx, ny = ny, ox = ox, oy = oy)
		) + theme(legend.position = "none")
    )
	}
    return(list(nx = nx, ny = ny, ox = ox, oy = oy))
}
