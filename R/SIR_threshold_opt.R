#' SIR optimally thresholded
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on (X,Y) with H slices, 
#' with thresholding of the matrix of interest by the *optimal* lambda parameter.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @param N.lambda The number of lambda to test. The N.lambda tested lambdas are 
#' between 0 and the maximum value of the interest matrix.
#' @param thresholding The thresholding method (choose between hard, soft)
#' @param graphic A boolean, set to TRUE to plot graphs 
#' @param output A boolean, set to TRUE to print informations
#' @return An object of class SIR.threshold.opt, with attributes:
#' \item{b.opt}{This is the optimal estimated EDR direction, which is the principal 
#' eigenvector of the interest matrix.}
#' \item{lambdas}{A vector that contains the tested lambdas}
#' \item{lambda.opt}{The optimal lambda}
#' \item{mat.b.th}{A matrix of size p*N.lambda that contains an estimation of beta 
#' in the columns for each lambda}
#' \item{N.lambda}{The number of lambda tested}
#' \item{vect.nb.zeros}{The number of 0 in the estimation of the vector beta}
#' \item{list.relevant.variables}{A list that contains the variables selected by 
#' the model}
#' @examples 
#' 
#' # Generate Data
#' set.seed(10)
#' n <- 200
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with hard thresholding
#' SIR_threshold_opt(Y,X,H=10,n_lambda=300,thresholding="hard",graphic=TRUE,output=TRUE)
#' 
#' # Apply SIR with soft thresholding
#' SIR_threshold_opt(Y,X,H=10,n_lambda=300,thresholding="soft",graphic=TRUE,output=TRUE)
#' @export
#' @importFrom strucchange breakpoints
SIR_threshold_opt <- function(Y, X, H = 10, n_lambda = 100, thresholding = "hard",
    graphic = TRUE, output = TRUE) {

    cl <- match.call()

    n <- nrow(X)
    p <- ncol(X)

    if (is.null(colnames(X))) {
        colnames(X) <- paste("X", 1:p, sep = "")
    }

    # Estimation de la direction des beta et de la matrice d'intérêt 
    # M1 = Sigma^-1 * Cov(Moyenne par tranche) avec la méthode SIR classique :
    res_SIR <- SIR(Y, X, H = 10, graphic = FALSE)
    b_sir <- res_SIR$b
    M1 <- res_SIR$M1

    # Création d'une liste de lambdas allant de 0 à la valeur absolue maximum 
    # de M1, avec un total de 100 valeurs.
    lambdas <- seq(0, max(abs(M1)), length.out = n_lambda + 1)[-(n_lambda + 1)]

    # Initialisation d'une matrice de 0 de taille N.lambda*p
    # pour reccueillir les estimations du vecteur b pour chacun des lambdas
    mat_b <- matrix(0, ncol = p, nrow = n_lambda)

    # Initialisation d'un vecteur de N.lambda NA pour receuillir le nombre de 
    # 0 après seuillage pour chaque valeur de lambda dans l'estimation du vecteur beta
    vect_nb_zeros <- rep(NA, n_lambda)

    # Initialisation d'un vecteur de N.lambda NA pour recueillir les cos carrés 
    # entre b_sir et b_seuillage pour chaque valeur de lambda
    vect_cosca <- rep(NA, n_lambda)

    # Initialisation de la liste des variables utiles
    list_relevant_variables <- list()

    # Application de la méthode SIR avec seuillage, avec les N.lambda valeurs
    # de lambdas
    for (i in 1:n_lambda) {

        # Récupération du résultat après application de la méthode SIR avec le lambda_i
        res_SIR_th <- SIR_threshold(Y, X, H = H, lambda = lambdas[i],
            thresholding = thresholding, graphic = FALSE)

        # Stockage de l'estimation du vecteur beta dans la ligne i de la matrice
        mat_b[i,] <- res_SIR_th$b

        # Stockage du nombre de 0 dans l'estimation du vecteur beta après seuillage pour
        # cette valeur de lambda
        vect_nb_zeros[i] <- res_SIR_th$nb_zeros

        # Stockage du cos carré entre b et b_seuillage pour cette valeur de lambda
        vect_cosca[i] <- res_SIR_th$cos_squared

        # Stocage des variables utiles pour cette valeur de lambda à l'indice i de 
        # la liste des variables utiles
        list_relevant_variables[[i]] <- res_SIR_th$list_relevant_variables
    }


    # Création d'un vecteur qui contient p valeurs : pour chaque variable x_i
    # est associé l'indice du lambda à partir duquel la variable devient inutile. 
    # Cet indice correspond donc également au nombre de lambda pour lesquels la variable
    # est utile.
    indices_useless_var <- colSums(mat_b / mat_b, na.rm = TRUE)
    names(indices_useless_var) = colnames(X)

    # On recherche ensuite un point de rupture dans la liste indice_useless_var ordonnée. 
    # Ce point de rupture fixe le nombre de variables à ne pas sélectionner, soit 
    # le nombre de variable où l'estimation de beta donne zero pour un lambda donné.
    fit_bp <- breakpoints(sort(indices_useless_var, decreasing = FALSE) ~ 1,
        breaks = 1, h = 2 / p)

    # Si le nombre de coeficients de l'estimation de beta égaux à 0 et qui sont égaux
    # au point de rupture, donc au nombre de variable à ne pas sélectionner : 
    if (length(which(vect_nb_zeros == fit_bp$breakpoints)) > 0) {
        # l'indice optimal de lambda est l'indice minimum où le nombre de 0 
        # apparaissant pour une valeur de lambda est égal au nombre de valeurs 
        # inutiles au point de rupture
        indice_opt <- min(which(vect_nb_zeros == fit_bp$breakpoints))
    } else {
        # Sinon, on somme le nombre de fois où le nombre de 0 qui apparait dans 
        # l'estimation de b pour chaque lambda est inférieur au nombre de valeurs 
        # inutiles au point de rupture + 1
        indice_opt <- sum(vect_nb_zeros < fit_bp$breakpoints) + 1
    }

    # On récupère le lambda qui correspond à cet indice
    lambda_opt <- lambdas[indice_opt]

    if (is.na(lambda_opt) == TRUE) {
        b <- SIR(Y, X, H = 10, graphic = FALSE)$b
        list_relevant_var <- colnames(X)
    }
    else {
        # Cas où on trouve bien un lambda optimal parmi les N.lambda
        # Si le nombre de variables inutiles est inférieur au nombre de variable total-1
        if (fit_bp$breakpoints < (p - 1)) {

            # Le beta optimal est la ligne de la matrice d'estimations de beta qui 
            # correspond à lambda.opt
            b <- mat_b[which(lambdas == lambda_opt),]
            # Conversion du beta en matrice à une ligne p colonnes
            b <- matrix(b, nrow = 1)
            # Renommage des colonnes
            colnames(b) <- colnames(X)
            # Récupération des variables utiles qui sont les colonnes du b.opt
            # où la valeur est différente de zero
            list_relevant_var <- colnames(b)[-which(b == 0)]


        }
        # Cas où la méthode Sparse SIR avec seuillage n'a pas permis de réaliser de la
        # sélection de variable
        else {
            b <- SIR(Y, X, H = 10, graphic = FALSE)$b
            list_relevant_var <- colnames(X)
        }
    }

    X_reduced <- X[, list_relevant_var, drop = FALSE]

    index_pred <- X %*% t(b)


    res <- list(b = b, lambdas = lambdas, lambda_opt = lambda_opt,
        mat_b = mat_b, n_lambda = n_lambda, vect_nb_zeros = vect_nb_zeros,
        fit_bp = fit_bp, indices_useless_var = indices_useless_var,
        vect_cosca = vect_cosca, list_relevant_variables = list_relevant_var,
        n = n, p = p, H = H, M1 = M1, thresholding = thresholding, call = cl,
        X_reduced = X_reduced, index_pred = index_pred, Y = Y)

    class(res) <- "SIR_threshold_opt"



    # Affichage de la sélection du lambda
    if (graphic == TRUE) {
        plot.SIR_threshold_opt.R(res)
    }

    return(res)
}
