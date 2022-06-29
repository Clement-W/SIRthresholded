#'  Vanilla SIR
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on (X,Y) with H slices.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @return An object of class SIR, with attributes:
#' \item{beta}{This is an estimated EDR direction, which is the principal 
#' eigenvector of the interest matrix.}
#' \item{M1}{The interest matrix.}
#' \item{eig.val}{The eigenvalues of the interest matrix.}
#' \item{eig.vect}{A matrix corresponding to the eigenvectors of the interest matrix.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
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
#' SIR(Y,X,H=10)
#' @export
SIR <- function(Y, X, H = 10, graphic = TRUE) {

    cl <- match.call()

    # Pour gérer le cas ou X ne contient qu'une variable
    if (is.null(dim(X)) || length(dim(X)) == 1) {
        X <- matrix(X, ncol = 1)
    }

    n <- nrow(X)
    p <- ncol(X)

    # Moyenne des échantillons pour chaque variable de X
    moy_X <- matrix(apply(X, 2, mean), ncol = 1)

    # Calcul de la matrice de covariance de X
    sigma <- var(X) * (n - 1) / n

    # Ordonne les données de X à partir des indices triés dans l'ordre croissant de Y
    X_ord <- X[order(Y),, drop = FALSE]

    # le vecteur nH contient le nombre d'observation pour chaque tranche H
    vect_nh <- rep(n %/% H, H) # %/% = integer division

    # si le nombre d'observations total est plus grand que le nb d'observation 
    # dans vect.nh pour prendre le reste de la division de n %/% H 
    if ((n - sum(vect_nh)) > 0) {
        # pour chaque observation non présente dans vect.nh
        for (i in 1:(n - sum(vect_nh))) {
            # on prend une tranche au hasard
            h <- sample(1:H, 1)
            # on ajoute 1 observation dans cette tranche
            vect_nh[h] <- vect_nh[h] + 1
        }
    }

    # vect.ph = proportion des yi tombant dans chaque tranche
    vect_ph <- vect_nh / n

    # séparation des n échantillons en H tranche
    # (on assigne un indice h à chaque échantillon)
    vect_h <- rep(1:H, times = vect_nh)

    # Matrice des moyennes par tranche des xi (matrice p*H)
    mat_mh <- matrix(0, ncol = H, nrow = p)

    # Pout chaque tranche
    for (h in 1:H) {
        # récupération des indices pour la tranche h
        index_slice_h <- which(vect_h == h)
        # récupération des x appartenant à la tranche h
        x_in_slice_h <- X_ord[index_slice_h,, drop = FALSE]
        # applique la moyenne aux échantillons de chaque feature de la tranche h
        mat_mh[, h] <- apply(x_in_slice_h, 2, mean)
    }

    # Estimation de la matrice de covariances des moyennes par tranche
    # (moy.X.etendue pour broadcaster l'opération - entre mat.mh et moy.X)
    moy_X_extended <- moy_X %*% matrix(rep(1, H), ncol = H)
    M <- (mat_mh - moy_X_extended) %*% diag(vect_ph) %*% t(mat_mh - moy_X_extended)

    # Matrice d'intérêt : inverse de la matrice de covariance de X multipliée par la
    # matrice de covariance des moyennes par tranche :(pxp) %*% (pxp)
    # matrice d'intérêt =  Sigma^-1 * Cov(Moyenne par tranche)
    sigma_inv <- solve(sigma)
    M1 <- sigma_inv %*% M

    # eigen(M1)$vectors donne une matrice pxp qui contient les p vecteurs propres de x 
    # dans chaque colonne. Les vecteurs propres sont retournés dans l'ordre décroissant
    # de leur valeur propre associée. Donc ici on récupère le vecteur propre associé 
    # à la plus grande valeur propre pour avoir b.est donc la première colonne
    res_SIR <- eigen(M1)
    eig_values <- Re(res_SIR$values)
    eig_vectors <- Re(res_SIR$vectors)
    b <- eig_vectors[, 1]

    # conversion en matrice à une ligne et p colonnes
    b <- matrix(b, nrow = 1)

    # Ajout des noms de colonnes
    if (!is.null(colnames(X))) {
        colnames(b) <- colnames(X)
    } else {
        colnames(b) <- paste("X", 1:p, sep = "")
    }
    # On a donc b.est l'estimation de la direction de Bet?a

    index_pred <- X %*% t(b)

    res <- list(b = b, M1 = M1, eig_val = eig_values,
        n = n, p = p, H = H, call = cl, index_pred = index_pred,
        Y = Y)
    class(res) <- "SIR"

    if (graphic) {
        plot.SIR(res)
    }

    return(res)
}