#' SIR thresholded
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on (X,Y) with H slices, with 
#' thresholding of the matrix of interest by the lambda parameter.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @param lambda The thresholding parameter
#' @param thresholding The thresholding method (choose between hard, soft)
#' @return An object of class SIR.threshold, with attributes:
#' \item{beta}{This is an estimated EDR direction, which is the principal eigenvector 
#' of the interest matrix.}
#' \item{M1_th}{The interest matrix thresholded.}
#' \item{eig.val}{The eigenvalues of the interest matrix thresholded.}
#' \item{eig.vect}{A matrix corresponding to the eigenvectors of the interest matrix.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' \item{nb.zeros}{The number of 0 in the estimation of the vector beta}
#' \item{list.relevant.variables}{A list that contains the variables selected by the
#' model}
#' \item{cos.squared}{The cosine squared between vanilla SIR and SIR thresholded}
#' @examples 
#' 
#' # Generate Data
#' set.seed(10)
#' n <- 500
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with hard thresholding
#' SIR_threshold(Y,X,H=10,lambda=0.2,thresholding="hard")
#' 
#' # Apply SIR with soft thresholding
#' SIR_threshold(Y,X,H=10,lambda=0.2,thresholding="soft")
#' @export
SIR_threshold <- function(Y, X, H = 10, lambda = 0, thresholding = "hard", graphic = TRUE) {

    cl <- match.call()

    n <- nrow(X)
    p <- ncol(X)

    if (is.null(colnames(X))) {
        colnames(X) <- paste("X", 1:p, sep = "")
    }

    # Estimation de la direction des beta et de la matrice d'intérêt 
    # Sigma^-1 * Cov(Moyenne par tranche) avec la méthode SIR classique :
    res_SIR <- SIR(Y, X, H = H, graphic = FALSE)
    b_sir <- res_SIR$b
    M1 <- res_SIR$M1

    # Seuillage de la matrice d'intérêt avec la méthode indiquée en paramètre
    if (thresholding == "soft") {
        M1_th <- do_soft_thresholding(M1, lambda = lambda)
    }
    if (thresholding == "hard") {
        M1_th <- do_hard_thresholding(M1, lambda = lambda)
    }

    res_eig <- eigen(M1_th)
    eig_values <- Re(res_eig$values)
    eig_vectors <- Re(res_eig$vectors)

    # Récupération du vecteur propre associé à la plus grande valeur propre
    # de la matrice d'intérêt seuillée
    b <- eig_vectors[, 1]

    # nombre de zéros présents dans l'estimation du b
    nb_zeros <- length(which(b == 0))

    # Calcul de la qualité de la corrélation entre les b estimés par la méthode classqiue
    # et par le sir thresholding
    cos_squared <- cosine_squared(b_sir, b)

    # Conversion en matrice en une ligne
    b <- matrix(b, nrow = 1)

    colnames(b) <- colnames(X)

    # Liste des variables utiles qui sont les colonnes de l'estimation de b
    # qui sont différentes de 0
    if (nb_zeros == 0) {
        list_relevant_variables <- colnames(X)
    } else {
        list_relevant_variables <- colnames(b)[-which(b == 0)]
    }

    X_reduced <- X[, list_relevant_variables, drop = FALSE]

    index_pred <- X %*% t(b)

    res <- list(b = b, M1 = M1_th, eig_val = eig_values, Y = Y,
        n = n, p = p, H = H, nb_zeros = nb_zeros, index_pred = index_pred,
        list_relevant_variables = list_relevant_variables, cos_squared = cos_squared,
        lambda = lambda, thresholding = thresholding, call = cl, X_reduced = X_reduced)
    class(res) <- "SIR_threshold"

    if (graphic) {
        plot.SIR_threshold(res)
    }

    return(res)

}