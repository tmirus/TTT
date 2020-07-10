#' Helper function for error handling
#' @details This function catches errors due to truncation of the lasso solution path and returns output for the best available lambda
#' @param lasso lasso model object you want to evaluate
#' @param lambda lambda value at which you are trying to evaluate
#' @return returns the nearest possible evaluation to your lambda
 
coef_error <- function(lasso,lambda){
    tmp <- try(genlasso::coef.genlasso(lasso, lambda=lambda),silent = TRUE)
    # if an error occured catch it, extract failed lambda from error message and recursively call 
    # coeff_error with slightly different lambda
    if (typeof(tmp)=='character'){
        new_lambda <- as.numeric(strsplit(strsplit(tmp[1],split='truncated at ')[[1]][2],split=',')[[1]][1])
        new_lambda <- new_lambda + 0.001
        warning(paste(lambda,' was not possible due to truncation, lambda used = ',new_lambda,sep = ''))
        tmp <- coef_error(lasso, lambda=new_lambda)
    }
    return(tmp)
}
