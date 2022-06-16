#'  Bootstrap SIR
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on B bootstraped sample of (X,Y) with H slices.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @param B The number of bootstrapped sample to draw
#' @return An object of class SIR.bootstrap, with attributes:
#' \item{beta}{This is an estimated EDR direction, which is the principal eigenvector of the interest matrix.}
#' \item{mat.b.est}{A matrix of size p*B that contains an estimation of beta in the columns for each bootstrapped sample}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' @examples 
#' # Generate Data
#' set.seed(10)
#' n <- 500
#' beta <- c(1,1,rep(0,8))
#' X <- rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply bootstrap SIR
#' SIR.bootstrap(Y,X,H=10,B=10)
#' @export
SIR.bootstrap <- function(Y, X, H = 10, B = 10) {
    p <- ncol(X)
    n <- nrow(X)

    # Initialisation de la matrice des estimations de beta_hat de taille pxB
    # avec l'estimation de beta_hat pour le nouvel échantillon après bootstrap dans
    # les colonnes de la matrice
    mat.b.est <- matrix(0, ncol = B, nrow = p)
    for (r in 1:B) {
        # indice des éléments à prendre dans X : on tire avec remise les indices des éléments de X
        # et on créé une liste qui contient ces indices (Bootstrap)
        indice <- sample(1:n, replace = TRUE)
        # on stock chaque estimation b dans les colonnes de la matrice d'estimation
        mat.b.est[, r] <- SIR(Y[indice], X[indice,], H = H)$beta
    }

    # on récupère le vecteur propre associé à la plus grande valeur propre
    # de la matrice d'estimations de beta_hat
    b.boot <- eigen(mat.b.est %*% t(mat.b.est), symmetric = TRUE)$vectors[, 1]

    # conversion en matrice à une ligne
    b.boot <- matrix(b.boot, nrow = 1)

    if (!is.null(colnames(X))) {
        colnames(b.boot) <- colnames(X)
    } else {
        colnames(b.boot) = paste("X", 1:p, sep = "")
    }

    res = list(beta = b.boot, mat.b.est = mat.b.est, n = n, p = p)
    class(res) = "SIR.bootstrap"

    return(res)
}