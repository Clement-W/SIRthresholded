#' SIR optimally thresholded on bootstraped replications
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on N bootstraped replications
#' of (X,Y) with H slices, with thresholding of the matrix of interest by the
#' *optimal* lambda parameter.
#' TODO: improve description
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @param N.lambda The number of lambda to test. The N.lambda tested lambdas are
#' between 0 and the maximum value of the interest matrix.
#' @param thresholding The thresholding method (choose between hard, soft)
#' @param Nb.replications The number of bootstraped replications of (X,Y) done to
#' estimate the model.
#' @param k Multiplication factor of the bootstrapped sample size
#' (default is 1 = keep the same size as original data)
#' @param graphic A boolean, set to TRUE to plot graphs
#' @param output A boolean, set to TRUE to print information
#' @return An object of class SIR.threshold.opt, with attributes:
#' \item{b.opt}{This is the optimal estimated EDR direction, which is the principal
#' eigenvector of the interest matrix.}
#' \item{lambdas}{A vector that contains the tested lambdas}
#' \item{lambda.opt}{The optimal lambda}
#' \item{Nb.var.selec.opt}{The optimal number of variables selected}
#' \item{list.relevant.variables}{A list that contains the variables selected
#' by the model}
#' @examples
#'
#' # Generate Data
#' set.seed(10)
#'  n <-  200
#' beta <- c(1,1,rep(0,18))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,20))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with hard thresholding
#' SIR.threshold.bootstrap(Y,X,H=10,N.lambda=300,thresholding="hard", Nb.replications=100,k=1,graphic=TRUE,output=TRUE)
#' 
#' # Apply SIR with soft thresholding
#' SIR.threshold.bootstrap(Y,X,H=10,N.lambda=300,thresholding="soft",Nb.replications=100,k=1,graphic=TRUE,output=TRUE)
#' @export
SIR.threshold.bootstrap <- function(Y, X, H = 10, thresholding = "hard",
    Nb.replications = 100, graphic = TRUE, output = TRUE, N.lambda = 100, k = 1) {

    cl <- match.call()

    # Sparse SIR avec N.lambda sur tout l'échantillon
    res.SparseSIR <- SIR.threshold.opt(Y, X, H = H, thresholding = thresholding,
        graph = FALSE, output = FALSE, N.lambda = N.lambda)

    p <- ncol(X)
    n <- nrow(X)

    if (is.null(colnames(X))) {
        colnames(X) = paste("X", 1:p, sep = "")
    }

    # Stockage du nombre de variable sélectionnée pour chaque réplication
    Nb.var.selec <- rep(NA, Nb.replications)
    # vecteur des noms de variables de X
    list.relevant.variables <- colnames(X)

    # Tableau permettant de stocker les variables sélectionnées pour chaque réplication
    liste <- list()

    # Vecteur qui contient le nombre de fois ou la variable à l'indice j a été retenue
    # dans une réplication
    effectif.var <- rep(0, p)

    # Pour chaque réplication
    for (replic in 1:Nb.replications) {
        # bootstrap 
        indice.boot <- sample(1:n, size = k * n, replace = TRUE)
        # On récupère les X et Y avec les indices pris par bootstrap
        X.boot <- X[indice.boot,]
        Y.boot <- Y[indice.boot]

        # Sparse.SIR sur l'echantillon créé par boostrap
        res.boot <- SIR.threshold.opt(Y.boot, X.boot, H = H, thresholding = thresholding,
            graph = FALSE, output = FALSE, N.lambda = N.lambda)
        # Stockage du nombre de variable sélectionnée (utiles) pour cette réplication
        Nb.var.selec[replic] <- length(res.boot$list.relevant.variables)

        # stockage des variables sélectionnées
        liste[[replic]] <- res.boot$list.relevant.variables

        if (output && replic %% 10 == 0) {
            print(paste("Replication n°", replic, "/", Nb.replications))
        }
    }

    # Pour chaque réplication
    for (i in 1:Nb.replications) {
        # Pour chaque variable
        for (j in 1:p) {
            # Si la variable à l'indice j a été retenue pour la ième réplication
            if ((sum(liste[[i]] == colnames(X)[j])) == 1) {
                # On incrémente le compteur d'occurrence de la variable j
                effectif.var[j] <- effectif.var[j] + 1
            }
        }
    }


    # Nombre de variable sélectionné optimale (le nombre de variable sélectionné 
    # qui est revenu le plus souvent parmi les réplications effectuées)
    Nb.var.selec.opt <-
        as.numeric(names(which(table(Nb.var.selec) == max(table(Nb.var.selec)))))

    # Nombre de zero optimal
    Nb.zeros.opt <- p - Nb.var.selec.opt

    # estimation du beta final en prenant le beta estimé sur tout l'échantillon par la
    # méthode SIR, au lambda à partir duquel le nombre de zero optimal apparaît
    b <- res.SparseSIR$mat.b.th[min(which(res.SparseSIR$vect.nb.zeros == Nb.zeros.opt)),]
    # Conversion du beta en matrice à une ligne p colonnes
    b <- matrix(b, nrow = 1)
    # Renommage des colonnes
    colnames(b) <- colnames(X)

    # Si on a bien réduit le nombre de variables
    if (Nb.zeros.opt > 0) {
        # mise à jour des variables utiles 
        list.relevant.variables <- colnames(X)[-which(b == 0)]
        # Si le nombre de zéros dans le b.opt final est à 0
        if (length(which(b == 0)) == 0) {
            # la liste de variables utile contient toutes les variables
            list.relevant.variables <- colnames(X)
        }
    }

    lambda.optim <-
        res.SparseSIR$lambdas[min(which(res.SparseSIR$vect.nb.zeros == Nb.zeros.opt))]

    res <- list(b = b, lambda.opt = lambda.optim,
        Nb.var.selec = Nb.var.selec, effectif.var = effectif.var, call = cl,
        Nb.var.selec.opt = Nb.var.selec.opt, list.relevant.variables =
        list.relevant.variables, n = n, p = p, H = H, Nb.replications =
        Nb.replications, thresholding = thresholding)

    class(res) <- "SIR.threshold.bootstrap"

    if (graphic == TRUE) {
        plot.SIR.threshold.bootstrap(res)
    }

    return(res)
}
