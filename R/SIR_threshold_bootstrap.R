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
#' SIR_threshold_bootstrap(Y,X,H=10,n_lambda=300,thresholding="hard", n_replications=100,k=1,graphic=TRUE,output=TRUE)
#' 
#' # Apply SIR with soft thresholding
#' SIR_threshold_bootstrap(Y,X,H=10,n_lambda=300,thresholding="soft",n_replications=100,k=1,graphic=TRUE,output=TRUE)
#' @export
SIR_threshold_bootstrap <- function(Y, X, H = 10, thresholding = "hard",
    n_replications = 50, graphic = TRUE, output = TRUE, n_lambda = 300, k = 2) {

    cl <- match.call()

    # Sparse SIR avec N.lambda sur tout l'échantillon
    res_SIR_th <- SIR_threshold_opt(Y, X, H = H, thresholding = thresholding,
        graph = FALSE, output = FALSE, n_lambda = n_lambda)

    p <- ncol(X)
    n <- nrow(X)

    if (is.null(colnames(X))) {
        colnames(X) <- paste("X", 1:p, sep = "")
    }

    # Stockage du nombre de variable sélectionnée pour chaque réplication
    nb_var_selec <- rep(NA, n_replications)
    # vecteur des noms de variables de X
    list_relevant_variables <- colnames(X)

    # Tableau permettant de stocker les variables sélectionnées pour chaque réplication
    liste <- list()

    # Vecteur qui contient le nombre de fois ou la variable à l'indice j a été retenue
    # dans une réplication
    effectif_var <- rep(0, p)

    mat_b <- matrix(0, ncol = p, nrow = n_replications)
    
    lambdas_opt_boot <- rep(0,n_replications)


    # Pour chaque réplication
    for (replic in 1:n_replications) {
        # bootstrap 
        indice_bootstrap <- sample(1:n, size = k * n, replace = TRUE)
        # On récupère les X et Y avec les indices pris par bootstrap
        X_boot <- X[indice_bootstrap,]
        Y_boot <- Y[indice_bootstrap]

        # Sparse.SIR sur l'echantillon créé par boostrap
        res_boot <- SIR_threshold_opt(Y_boot, X_boot, H = H, thresholding = thresholding,
            graph = FALSE, output = FALSE, n_lambda = n_lambda)
        # Stockage du nombre de variable sélectionnée (utiles) pour cette réplication
        nb_var_selec[replic] <- length(res_boot$list_relevant_variables)

        # stockage des variables sélectionnées
        liste[[replic]] <- res_boot$list_relevant_variables

        mat_b[replic,] <- res_boot$b
        
        lambdas_opt_boot[replic] <- res_boot$lambda_opt

        if (output && replic %% 5 == 0) {
            print(paste("Replication n°", replic, "/", n_replications))
        }
    }

    # Pour chaque réplication
    for (i in 1:n_replications) {
        # Pour chaque variable
        for (j in 1:p) {
            # Si la variable à l'indice j a été retenue pour la ième réplication
            if ((sum(liste[[i]] == colnames(X)[j])) == 1) {
                # On incrémente le compteur d'occurrence de la variable j
                effectif_var[j] <- effectif_var[j] + 1
            }
        }
    }


    # Nombre de variable sélectionné optimale (le nombre de variable sélectionné 
    # qui est revenu le plus souvent parmi les réplications effectuées)
    nb_var_selec_opt <-
        as.numeric(names(which(table(nb_var_selec) == max(table(nb_var_selec)))))[1]

    # Nombre de zero optimal
    nb_zeros_opt <- p - nb_var_selec_opt

    # estimation du beta final en prenant le beta estimé sur tout l'échantillon par la
    # méthode SIR, au lambda à partir duquel le nombre de zero optimal apparaît
    b <- res_SIR_th$mat_b[min(which(res_SIR_th$vect_nb_zeros == nb_zeros_opt)),]
    # Conversion du beta en matrice à une ligne p colonnes
    b <- matrix(b, nrow = 1)
    # Renommage des colonnes
    colnames(b) <- colnames(X)

    # Si on a bien réduit le nombre de variables
    if (nb_zeros_opt > 0) {
        # mise à jour des variables utiles 
        list_relevant_variables <- colnames(X)[-which(b == 0)]
        # Si le nombre de zéros dans le b.opt final est à 0
        if (length(which(b == 0)) == 0) {
            # la liste de variables utile contient toutes les variables
            list_relevant_variables <- colnames(X)
        }
    }

    X_reduced <- X[, list_relevant_variables, drop = FALSE]

    lambda_optim <-
        res_SIR_th$lambdas[min(which(res_SIR_th$vect_nb_zeros == nb_zeros_opt))]

    res <- list(b = b, lambda_opt = lambda_optim,
        nb_var_selec = nb_var_selec, effectif_var = effectif_var, call = cl,
        nb_var_selec_opt = nb_var_selec_opt, list_relevant_variables =
        list_relevant_variables, n = n, p = p, H = H, n_replications =
        n_replications, thresholding = thresholding, X_reduced = X_reduced, mat_b = mat_b,
        lambdas_opt_boot=lambdas_opt_boot)

    class(res) <- "SIR_threshold_bootstrap"

    if (graphic == TRUE) {
        plot.SIR_threshold_bootstrap(res)
    }

    return(res)
}
