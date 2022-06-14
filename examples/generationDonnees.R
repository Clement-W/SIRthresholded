library(mvtnorm)
library(strucchange)

#========== Génération de données 1 ========= 
# X ~ loi normale multivariée
# Epsilon ~ Loi normale de variance 10
# Fonction de lien  g(x)=x^3
# Modele : y = g(Bx) + Epsilon

# n = nombre d'observations
# p = nombre de dimensions
# p1 = nombre de coefficients non nuls dans Beta (variables utiles)
generationDonnees1 <- function(n, p, p1 = 5, graphic = FALSE) {

    # Données X complétement décorrélées
    Sigma <- diag(p)

    # multivariate normal distribution X  avec n observations et une matrice de covariance Sigma
    X <- rmvnorm(n, sigma = Sigma)
    dimnames(X) <- list(1:n, paste("X", 1:p, sep = ""))

    # Erreur epsilon de loi normale et de variance 10
    eps <- rnorm(n, sd = 10)

    # Création du vecteur beta
    beta.simple <- rep(1, p1) # liste de taille p1
    beta.simple <- matrix(beta.simple, ncol = 1) # matrice colonne de taille p1
    beta <- c(beta.simple, rep(0, p - p1)) # on ajoute des 0 à la suite pour avoir la même taille que X
    beta <- matrix(beta, ncol = 1) # on retransforme le beta en matrice colonne

    # renommage des lignes de beta
    rownames(beta) <- colnames(X)
    rownames(beta.simple) <- colnames(X)[1:p1]

    # Vecteur y = g(X*Beta)+eps avec g la fonction cubique
    Y <- ((X %*% beta) ** 3) + eps

    # Ratio des variances pour montrer à quel point le Y est bruité
    ratio <- var(eps) / var(Y)

    if (graphic == TRUE) {
        plot(X %*% beta, Y, xlab = "true index")
        title(paste("ratio V(epsilon)/V(y)=", round(ratio, digits = 3)))
    }

    return(list(Y = Y, X = X, beta = beta, beta.simple = beta.simple, ratio = ratio))
}


#========== Génération de données 2 ==========
# X ~ loi normale multivariée
# Epsilon ~ Loi normale de variance dépendente du ratio passé en paramètre
# Fonction de lien  g(x)=x^3
# Modele : y = g(Bx) + Epsilon


# n = nombre d'observations
# p = nombre de dimensions
# p1 = nombre de coefficients non nuls dans Beta (variables utiles)
# RATIO = le ratio souhaité entre la variance de l'erreur epsilon et la variance de Y (var(eps)/var(Y))
# SigmaIdentite = TRUE Si la covariance de x est décrollée, FALSE sinon
# fct = la fonction de lien g utilisée
generationDonneesAvecRatio <- function(n, p, p1 = 5, graphic = FALSE, RATIO = 0.05, correlationX = FALSE, fct = "cubique") {


    if (!correlationX) {
        # Données X complétement décorrélées
        Sigma <- diag(p)
    } else {
        # Corrélation auto-regressive
        rho = 0.1
        Sigma <- matrix(0, p, p)
        for (j in 1:p) {
            for (k in 1:p) {
                Sigma[j, k] <- rho ^ (abs(k - j))
            }
        }
    }

    # Distribution normale multivarée X avec n observations et une matrice de covariance Sigma
    X <- rmvnorm(n, sigma = Sigma)
    dimnames(X) <- list(1:n, paste("X", 1:p, sep = ""))

    # Création du vecteur beta
    beta.simple <- rep(1, p1) # liste de taille p1
    #beta.simple[1:5] = 0.5

    #beta.simple <- (1:p1)/p1 
    beta.simple <- matrix(beta.simple, ncol = 1) # matrice colonne de taille p1
    beta <- c(beta.simple, rep(0, p - p1)) # on ajoute des 0 à la suite pour avoir la même taille que X
    beta <- matrix(beta, ncol = 1) # on retransforme le beta en matrice colonne

    # renommage des lignes de beta
    rownames(beta) <- colnames(X)
    rownames(beta.simple) <- colnames(X)[1:p1]

    # Création de Y, avec le ratio (var(eps)/var(Y)) supérieur
    # ou égal à l'argument RATIO passé en parmaètre :

    # Tant que le ratio de variance entre epsilon est inférieur ou égal à 
    # l'argument RATIO et que moins de 1000 itérations ont été effectuées,
    # on incrémente l'écart type pour atteindre le RATIO attendu

    ecart.type <- 0 # écart type initial
    cpt <- 1
    ratio <- -Inf
    while (cpt < 1000 && ratio < RATIO) {
        ecart.type <- ecart.type + 0.1
        cpt <- cpt + 1

        # Création du vecteur d'erreur epsilon de loi normale et de variance ecart.type
        eps <- rnorm(n, sd = ecart.type * sqrt(sum(Sigma ^ 2) ^ 1.5))

        if (fct == "cubique") {
            Y <- fct.cubique(X, beta, eps)
        }
        else if (fct == "exp") {
            Y <- fct.exp(X, beta, eps)
        }
        else if (fct == "exp-heteroscedastique") {
            Y <- fct.exp.heteroscedastique(X, beta, eps)
        }
        else if (fct == "lineaire") {
            Y <- fct.lineaire(X, beta, eps)
        }
        else if (fct == "sinexp") {
            Y <- fct.sinexp(X, beta, eps)
        }
        else {
            stop(paste("Error, function", fct, " is not defined"))
        }

        # Calcul du ratio de variance entre epsilon et Y
        ratio <- var(eps) / var(Y)
    }

    if (graphic == TRUE) {
        plot(X %*% beta, Y, xlab = "true index")
        title(paste("ratio V(epsilon)/V(y)=", round(ratio, digits = 3)))
    }

    return(list(Y = Y, X = X, beta = beta, beta.simple = beta.simple, ratio = ratio, cpt = cpt))
}

# Retourne un vecteur y = g(X*Beta)+eps avec g la fonction cubique
fct.cubique <- function(X, beta, eps) {
    return(((X %*% beta) ** 3) + eps)
}

# Retourne un vecteur y = g(X*Beta)+eps avec g la fonction exponentielle
fct.exp <- function(X, beta, eps) {
    return((exp(X %*% beta)) + eps)
}

# Retourne un vecteur y = g(X*Beta+eps) avec g la fonction exponentielle qui prend également l'erreur
fct.exp.heteroscedastique <- function(X, beta, eps) {
    return((exp(X %*% beta) + eps))
}

# Retourne un vecteur y = g(X*Beta+eps) avec g la fonction linéaire x*beta + eps
fct.lineaire <- function(X, beta, eps) {
    return((X %*% beta) + eps)
}

# Retourne un vecteur y = g(X*Beta+eps) avec g la sinus * exponentielle
fct.sinexp <- function(X, beta, eps) {
    return(sin(X %*% beta) * exp(X %*% beta) + eps)
}