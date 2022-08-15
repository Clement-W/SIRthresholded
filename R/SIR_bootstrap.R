#'  Bootstrap SIR
#'
#' Apply a single-index \eqn{SIR} on \eqn{B} bootstraped samples of \eqn{(X,Y)} with \eqn{H} slices. 
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices (default is 10).
#' @param B The number of bootstrapped samples to draw (default is 10).
#' @param graph A boolean that must be set to true to display graphics (default is TRUE).
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "eigvals" Plot the eigen values of the matrix of interest.
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y.
#'   \item "" Plot every graphs (default).
#' }
#' @return An object of class SIR_bootstrap, with attributes:
#' \item{b}{This is an estimated EDR direction, which is the principal eigenvector
#' of the interest matrix.}
#' \item{mat_b}{A matrix of size p*B that contains an estimation of beta in 
#' the columns for each bootstrapped sample.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' \item{H}{The chosen number of slices.}
#' \item{call}{Unevaluated call to the function.}
#' \item{index_pred}{The index b'X estimated by SIR.}
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
#' # Apply bootstrap SIR
#' SIR_bootstrap(Y, X, H = 10, B = 10)
#' @export
SIR_bootstrap <- function(Y, X, H = 10, B = 10, graph = TRUE, choice = "") {

    cl <- match.call()

    p <- ncol(X)
    n <- nrow(X)

    # Initialization of a matrix of size pxB  that will contain the estimation 
    # of b for each the new sample after bootstrapping in the columns of the matrix
    mat_b <- matrix(0, ncol = B, nrow = p)
    for (r in 1:B) {
        # index of the samples to be taken in X : we draw with discount the indices
        # of X samples and store it in a vector (Bootstrap)
        indice <- sample(1:n, replace = TRUE)
        # store each estimate b in the columns of mat_b
        mat_b[, r] <- SIR(Y[indice], X[indice,], H = H, graph = FALSE)$b
    }

    # we recover the eigenvector associated to the largest eigenvalue of mat_b
    e <- eigen(mat_b %*% t(mat_b), symmetric = TRUE)
    b <- e$vectors[, 1]
    eig_values <- Re(e$values)

    # convert it into a one line matrix
    b <- matrix(b, nrow = 1)

    if (!is.null(colnames(X))) {
        colnames(b) <- colnames(X)
    } else {
        colnames(b) <- paste("X", 1:p, sep = "")
    }

    # compute the estimated index
    index_pred <- X %*% t(b)

    res <- list(b = b, mat_b = mat_b, n = n, p = p, H = H, call = cl,
    index_pred = index_pred, eig_values = eig_values, Y = Y)
    class(res) <- "SIR_bootstrap"

    if (graph) {
        plot.SIR_bootstrap(res, choice = choice)
    }

    return(res)
}