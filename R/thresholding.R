#'  Soft thresholding
#'
#' Apply a soft thresholding of parameter lambda to the matrix M
#' @param M The matrix that will be thresholded
#' @param lambda The parameter of the tresholding
#' @return The matrix M thresholded
#' @examples 
#' M = matrix(rnorm(36),nrow=6)
#' M1 = do_soft_thresholding(M = M, lambda = 0.4)
do_soft_thresholding <- function(M, lambda) {
    abs_m <- abs(M)
    M_th <- abs_m - lambda
    M_th[which(M_th < 0)] <- 0
    return(M_th * sign(M))
}

#'  Hard thresholding
#'
#' Apply a hard thresholding of parameter lambda to the matrix M
#' @param M The matrix that will be thresholded
#' @param lambda The parameter of the tresholding
#' @return The matrix M thresholded
#' @examples 
#' M = matrix(rnorm(36),nrow=6)
#' M1 = do_hard_thresholding(M = M, lambda = 0.4)
do_hard_thresholding <- function(M, lambda) {
    M.hard <- M
    abs_m <- abs(M)
    M_th <- abs_m - lambda
    M.hard[which(M_th < 0)] <- 0
    return(M.hard)
}

#TODO: Add SCAD thresholding