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

#TODO: Add SCAD thresholding

do_scad_thresholding <- function(M,lambda,a) {
    M.scad <- M
    abs_m <- abs(M)
    index_cond1 = which(abs_m <= 2*lambda)
    M.scad[index_cond1] = sign(M[index_cond1]) * ifelse((abs_m[index_cond1] - lambda) > 0,abs_m[index_cond1] - lambda,0) 
    index_cond2 = which((abs_m > 2*lambda) & (abs_m <= a*lambda))
    M.scad[index_cond2] = -sign(M[index_cond2]) * a * lambda + (a - 1) * M[index_cond2]
    return(M.scad)
}
