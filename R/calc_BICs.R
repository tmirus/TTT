#' Returns useful information about a fusedlasso model fit
#' This function was adapted from Martin Fahrenberger's sflST package 
#' (https://github.com/Martin-Fahrenberger/sflST)
#' 
#' @details This function calculates the RSS, BIC and explained variance of a genlasso model
#' @param full_model takes an object of the class fusedlasso
#' @return returns a vector consisting of the df, lambda, rss, BIC and ratio of explained variation dor each lambda.

calc_BICs <- function(full_model){
    summary_lassos <- genlasso::summary.genlasso(full_model)
    N_BIC <- length(full_model$y)
    SDS_BIC <- var(full_model$y)
    ev_ratio <- summary_lassos[,3]/(N_BIC*SDS_BIC)
    BIC <- ev_ratio + log(N_BIC)/N_BIC*summary_lassos[,1]
    
    return(cbind(summary_lassos,BIC,ev_ratio))
}
