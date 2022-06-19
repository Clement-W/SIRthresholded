#'  Vanilla SIR
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on (X,Y) with H slices.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @return An object of class SIR, with attributes:
#' \item{beta}{This is an estimated EDR direction, which is the principal eigenvector of the interest matrix.}
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
SIR <- function(Y, X, H = 10) {

    # Pour gérer le cas ou X ne contient qu'une variable
    if (is.null(dim(X)) || length(dim(X)) == 1) {
        X = matrix(X, ncol = 1)
    }

    n <- nrow(X)
    p <- ncol(X)

    # Moyenne des échantillons pour chaque variable de X
    moy.X <- matrix(apply(X, 2, mean), ncol = 1)

    # Calcul de la matrice de covariance de X
    Sigma <- var(X) * (n - 1) / n

    # Ordonne les données de X à partir des indices triés dans l'ordre croissant de Y
    X.ord <- X[order(Y),, drop = FALSE]

    # le vecteur nH contient le nombre d'observation pour chaque tranche H
    vect.nh <- rep(n %/% H, H) # %/% = integer division

    # si le nombre d'observations total est plus grand que le nb d'observation dans vect.nh
    # pour prendre le reste de la division de n %/% H 
    if ((n - sum(vect.nh)) > 0) {
        # pour chaque observation non présente dans vect.nh
        for (i in 1:(n - sum(vect.nh))) {
            # on prend une tranche au hasard
            h <- sample(1:H, 1)
            # on ajoute 1 observation dans cette tranche
            vect.nh[h] <- vect.nh[h] + 1
        }
    }

    # vect.ph = proportion des yi tombant dans chaque tranche
    vect.ph <- vect.nh / n

    # séparation des n échantillons en H tranche
    # (on assigne un indice h à chaque échantillon)
    vect.h <- rep(1:H, times = vect.nh)

    # Matrice des moyennes par tranche des xi (matrice p*H)
    mat.mh <- matrix(0, ncol = H, nrow = p)

    # Pout chaque tranche
    for (h in 1:H) {
        # récupération des indices pour la tranche h
        indicesTrancheH = which(vect.h == h)
        # récupération des x appartenant à la tranche h
        xDansTrancheH = X.ord[indicesTrancheH,, drop = FALSE]
        # applique la moyenne aux échantillons de chaque feature de la tranche h
        mat.mh[, h] <- apply(xDansTrancheH, 2, mean)
    }

    # Estimation de la matrice de covariances des moyennes par tranche
    moy.X.etendue = moy.X %*% matrix(rep(1, H), ncol = H) # Pour broadcaster l'opération - entre mat.mh et moy.X
    M <- (mat.mh - moy.X.etendue) %*% diag(vect.ph) %*% t(mat.mh - moy.X.etendue)

    # Matrice d'intérêt : inverse de la matrice de covariance de X multipliée par la
    # matrice de covariance des moyennes par tranche :(pxp) %*% (pxp)
    #  matrice d'intérêt =  Sigma^-1 * Cov(Moyenne par tranche)
    SigmaInv <- solve(Sigma)
    M1 <- SigmaInv %*% M

    # eigen(M1)$vectors donne une matrice pxp qui contient les p vecteurs propres de x dans chaque colonne
    # Les vecteurs propres sont retournés dans l'ordre décroissant de leur valeur propre associée.
    # Donc ici on récupère le vecteur propre associé à la plus grande valeur propre pour avoir b.est
    # donc la première colonne
    rSIR = eigen(M1)
    eig.values = Re(rSIR$values)
    eig.vectors = Re(rSIR$vectors)
    b.est <- eig.vectors[, 1]

    # conversion en matrice à une ligne et p colonnes
    b.est <- matrix(b.est, nrow = 1)

    # Ajout des noms de colonnes
    if (!is.null(colnames(X))) {
        colnames(b.est) <- colnames(X)
    } else {
        colnames(b.est) = paste("X", 1:p, sep = "")
    }
    # On a donc b.est l'estimation de la direction de Beta

    res = list(beta = b.est, M1 = M1, eig.val = eig.values, eig.vect = eig.vectors, n = n, p = p,H=H)
    class(res) = "SIR"

    return(res)
}

