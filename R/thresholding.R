# Apply a soft thresholding of parameter lambda to the matrix M
do_soft_thresholding <- function(M, lambda) {
    abs_m <- abs(M)
    M_th <- abs_m - lambda
    M_th[which(M_th < 0)] <- 0
    return(M_th * sign(M))
}

# Apply a hard thresholding of parameter lambda to the matrix M
do_hard_thresholding <- function(M, lambda) {
    M.hard <- M
    abs_m <- abs(M)
    M_th <- abs_m - lambda
    M.hard[which(M_th < 0)] <- 0
    return(M.hard)
}
