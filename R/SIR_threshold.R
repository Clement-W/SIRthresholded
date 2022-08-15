#' SIR threshold
#'
#' Apply a single-index \eqn{SIR} on \eqn{(X,Y)} with \eqn{H} slices, with a parameter \eqn{\lambda} which 
#' apply a soft/hard thresholding to the interest matrix \eqn{\widehat{\Sigma}_n^{-1}\widehat{\Gamma}_n}.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices (default is 10).
#' @param lambda The thresholding parameter (default is 0).
#' @param thresholding The thresholding method to choose between hard and soft (default is hard).
#' @param graph A boolean that must be set to true to display graphics (default is TRUE).
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "eigvals" Plot the eigen values of the matrix of interest.
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y.
#'   \item "" Plot every graphs (default).
#' }
#' @return An object of class SIR_threshold, with attributes:
#' \item{b}{This is an estimated EDR direction, which is the principal eigenvector 
#' of the interest matrix.}
#' \item{M1}{The interest matrix thresholded.}
#' \item{eig_val}{The eigenvalues of the interest matrix thresholded.}
#' \item{eig_vect}{A matrix corresponding to the eigenvectors of the interest matrix.}
#' \item{Y}{The response vector.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' \item{H}{The chosen number of slices.}
#' \item{nb.zeros}{The number of 0 in the estimation of the vector beta.}
#' \item{index_pred}{The index Xb' estimated by SIR.}
#' \item{list.relevant.variables}{A list that contains the variables selected by the
#' model.}
#' \item{cos_squared}{The cosine squared between vanilla SIR and SIR thresholded.}
#' \item{lambda}{The thresholding parameter used.}
#' \item{thresholding}{The thresholding method used.}
#' \item{call}{Unevaluated call to the function.}
#' \item{X_reduced}{The X data restricted to the variables selected by the model.
#' It can be used to estimate a new SIR model on the relevant variables to improve
#' the estimation of b.}
#' @examples 
#' # Generate Data
#' set.seed(10)
#' n <- 500
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with hard thresholding
#' SIR_threshold(Y, X, H = 10, lambda = 0.2, thresholding = "hard")
#' @export
SIR_threshold <- function(Y, X, H = 10, lambda = 0, thresholding = "hard", graph = TRUE, choice = "") {

    cl <- match.call()

    # Ensure that X and Y are matrices
    X = ensure_matrix(X)
    Y = ensure_matrix(Y)

    n <- nrow(X)
    p <- ncol(X)

    if (is.null(colnames(X))) {
        colnames(X) <- paste("X", 1:p, sep = "")
    }

    # Estimation of b and the interest matrix with the classic SIR method
    res_SIR <- SIR(Y, X, H = H, graph = FALSE)
    b_sir <- res_SIR$b
    M1 <- res_SIR$M1

    # Thresholding of the interest matrix
    if (thresholding == "soft") {
        M1_th <- do_soft_thresholding(M1, lambda = lambda)
    }
    if (thresholding == "hard") {
        M1_th <- do_hard_thresholding(M1, lambda = lambda)
    }

    # Compute eigenvalues and eigenvectors of the thresholded interest matrix
    res_eig <- eigen(M1_th)
    eig_values <- Re(res_eig$values)
    eig_vectors <- Re(res_eig$vectors)


    # Get the eigenvector associated to the greatest eigevalue of the 
    # thresholded interest matrix
    b <- eig_vectors[, 1]

    # Number of zeros in b (number of useless variables)
    nb_zeros <- length(which(b == 0))

    # Compute the quality of the correlation between b_sir and b_thresholded_sir
    cos_squared <- cosine_squared(b_sir, b)

    # Convert b into a one-line matrix
    b <- matrix(b, nrow = 1)
    colnames(b) <- colnames(X)

    # The list of relevant variables are the columns of b that are different to 0
    if (nb_zeros == 0) {
        list_relevant_variables <- colnames(X)
    } else {
        list_relevant_variables <- colnames(b)[-which(b == 0)]
    }

    # Create the X reduced variable by restricting X to the relevant variables.
    X_reduced <- X[, list_relevant_variables, drop = FALSE]

    index_pred <- X %*% t(b)

    res <- list(b = b, M1 = M1_th, eig_val = eig_values, Y = Y,
        n = n, p = p, H = H, nb_zeros = nb_zeros, index_pred = index_pred,
        list_relevant_variables = list_relevant_variables, cos_squared = cos_squared,
        lambda = lambda, thresholding = thresholding, call = cl, X_reduced = X_reduced)
    class(res) <- "SIR_threshold"

    if (graph) {
        plot.SIR_threshold(res, choice = choice)
    }

    return(res)

}