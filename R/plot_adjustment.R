#' @export
plot_adjustment <- function(ids, img, nx = 35, ny = 33, ox = 0, oy = 0){
    plot(spatial_plot(rownames(ids), ids, rep(1, nrow(ids)), img, "discrete", nx, ny, ox, oy))
    return(list(nx = nx, ny = ny, ox = ox, oy = oy))
}