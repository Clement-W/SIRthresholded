#'  Classic SIR
#'
#' Apply a single-index \eqn{SIR} on \eqn{(X,Y)} with \eqn{H} slices. This function allows to obtain an 
#' estimate of a basis of the \eqn{EDR} (Effective Dimension Reduction) space via the eigenvector 
#' \eqn{\hat{b}} associated with the largest nonzero eigenvalue of the matrix of interest
#' \eqn{\widehat{\Sigma}_n^{-1}\widehat{\Gamma}_n}. Thus, \eqn{\hat{b}} is an \eqn{EDR} direction.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices (default is 10).
#' @param graph A boolean that must be set to true to display graphics (default is TRUE).
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "eigvals" Plot the eigen values of the matrix of interest.
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y.
#'   \item "" Plot every graphs. (default)
#' }
#' @return An object of class SIR, with attributes:
#' \item{b}{This is an estimated EDR direction, which is the principal 
#' eigenvector of the interest matrix.}
#' \item{M1}{The interest matrix.}
#' \item{eig_val}{The eigenvalues of the interest matrix.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' \item{H}{The chosen number of slices.}
#' \item{call}{Unevaluated call to the function.}
#' \item{index_pred}{The index Xb' estimated by SIR.}
#' \item{Y}{The response vector.}
#' @examples
#' # Generate Data
#' set.seed(10)
#' n <- 500
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR
#' SIR(Y, X, H = 10)
#' @export
#' @importFrom stats var
SIR <- function(Y, X, H = 10, graph = TRUE, choice = "") {

    cl <- match.call()

    # Ensure that X and Y are matrices
    X = ensure_matrix(X)
    Y = ensure_matrix(Y)

    n <- nrow(X)
    p <- ncol(X)

    # Mean of each variables of the samples X
    moy_X <- matrix(apply(X, 2, mean), ncol = 1)

    # Compute the covariance matrix of X
    sigma <- var(X) * (n - 1) / n

    # Order the data of X from the indices sorted in ascending order of Y
    X_ord <- X[order(Y),, drop = FALSE]

    # the vector nH contains the number of observations for each slice H
    vect_nh <- rep(n %/% H, H) # %/% = integer division

    # Take into account the rest of the integer division:
    rest = n - sum(vect_nh)
    if (rest > 0) {
        # For each sample not in vect_nh
        h = 1
        for (i in 1:rest) {
            # add a sample in the slice h
            vect_nh[h] <- vect_nh[h] + 1

            # Increment h (which will be < H because rest < H)
            h <- h + 1
        }
    }

    # vect_ph = ratio of y_i falling in each slice
    vect_ph <- vect_nh / n


    # Separate the n samples in H slices by assigning an index h to each sample
    vect_h <- rep(1:H, times = vect_nh)
    # vect_h = 1 1 1 1 1 .... 1 2 2 2 2 ....... 10 10 10 10 10 10 if H=10

    #  matrix of size p*H that will contains the mean of the x_i by slice
    mat_mh <- matrix(0, ncol = H, nrow = p)

    # For each slice, compute the mean of the x_i
    for (h in 1:H) {
        # Get the indices in the slice h
        index_slice_h <- which(vect_h == h)
        # Get the x belonging to the slice h
        x_in_slice_h <- X_ord[index_slice_h,, drop = FALSE]
        # Compute the mean for each feature of x in the slice h
        mat_mh[, h] <- apply(x_in_slice_h, 2, mean)
    }

    # Estimation of the covariance matrix of the means per slice
    # (moy_X_extended is created to broadcast the operation - beween mat_mh and moy_X)
    moy_X_extended <- moy_X %*% matrix(rep(1, H), ncol = H)
    M <- (mat_mh - moy_X_extended) %*% diag(vect_ph) %*% t(mat_mh - moy_X_extended)


    # Interest matrix : inverse of the covariance matrix of X multiplied by the
    # covariance matrix of slice means : (pxp) %*% (pxp)
    # Interest matrix = Sigma^-1 * Cov(Slice mean)
    sigma_inv <- solve(sigma)
    M1 <- sigma_inv %*% M

    # eigen(M1)$vectors gives a matrix pxp which contains the p eigenvectors of x 
    # in each column. The eigenvectors are returned in the decreasing order
    # of their associated eigenvalue. So here we get the eigenvector associated with 
    # to the largest eigenvalue to have b thus the first column

    res_SIR <- eigen(M1)
    eig_values <- Re(res_SIR$values)
    eig_vectors <- Re(res_SIR$vectors)
    b <- eig_vectors[, 1]

    # convert into a one line and p columns matrix
    b <- matrix(b, nrow = 1)

    # Add column names 
    if (!is.null(colnames(X))) {
        colnames(b) <- colnames(X)
    } else {
        colnames(b) <- paste("X", 1:p, sep = "")
    }

    # Compute the index Xb'
    index_pred <- X %*% t(b)

    res <- list(b = b, M1 = M1, eig_val = eig_values,
        n = n, p = p, H = H, call = cl, index_pred = index_pred,
        Y = Y)
    class(res) <- "SIR"

    if (graph) {
        plot.SIR(res, choice = choice)
    }

    return(res)
}