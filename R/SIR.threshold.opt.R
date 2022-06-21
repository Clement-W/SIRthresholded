#' SIR optimally thresholded
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on (X,Y) with H slices, with thresholding of the matrix of interest by the *optimal* lambda parameter.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @param N.lambda The number of lambda to test. The N.lambda tested lambdas are between 0 and the maximum value of the interest matrix.
#' @param thresholding The thresholding method (choose between hard, soft)
#' @param graphic A boolean, set to TRUE to plot graphs 
#' @param output A boolean, set to TRUE to print informations
#' @return An object of class SIR.threshold.opt, with attributes:
#' \item{b.opt}{This is the optimal estimated EDR direction, which is the principal eigenvector of the interest matrix.}
#' \item{lambdas}{A vector that contains the tested lambdas}
#' \item{lambda.opt}{The optimal lambda}
#' \item{mat.b.th}{A matrix of size p*N.lambda that contains an estimation of beta in the columns for each lambda}
#' \item{N.lambda}{The number of lambda tested}
#' \item{vect.nb.zeros}{The number of 0 in the estimation of the vector beta}
#' \item{list.relevant.variables}{A list that contains the variables selected by the model}
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
#' SIR.threshold.opt(Y,X,H=10,N.lambda=300,thresholding="hard",graphic=TRUE,output=TRUE)
#' 
#' # Apply SIR with soft thresholding
#' SIR.threshold.opt(Y,X,H=10,N.lambda=300,thresholding="soft",graphic=TRUE,output=TRUE)
#' @export
SIR.threshold.opt <- function(Y, X, H = 10, N.lambda = 100, thresholding = "hard", graphic = TRUE, output = TRUE) {

    n <- nrow(X)
    p <- ncol(X)

    if (is.null(colnames(X))) {
        colnames(X) = paste("X", 1:p, sep = "")
    }

    # Estimation de la direction des beta et de la matrice d'intérêt  M1 = Sigma^-1 * Cov(Moyenne par tranche)
    # avec la méthode SIR classique :
    resSIR <- SIR(Y, X, H = 10)
    b.clas <- resSIR$beta
    M1 <- resSIR$M1

    # Création d'une liste de lambdas allant de 0 à la valeur absolue maximum de M1, avec un total
    # de 100 valeurs.
    lambdas <- seq(0, max(abs(M1)), length.out = N.lambda + 1)[-(N.lambda + 1)]

    # Initialisation d'une matrice de 0 de taille N.lambda*p
    # pour reccueillir les estimations du vecteur b pour chacun des lambdas
    mat.b.th <- matrix(0, ncol = p, nrow = N.lambda)

    # Initialisation d'un vecteur de N.lambda NA pour receuillir le nombre de 0 après seuillage
    # pour chaque valeur de lambda dans l'estimation du vecteur beta
    vect.nb.zeros <- rep(NA, N.lambda)

    # Initialisation d'un vecteur de N.lambda NA pour recueillir les cos carrés entre b_sir et b_seuillage
    # pour chaque valeur de lambda
    vect.cosca <- rep(NA, N.lambda)

    # Initialisation de la liste des variables utiles
    list.relevant.variables <- list()

    # Application de la méthode SIR avec seuillage, avec les N.lambda valeurs
    # de lambdas
    for (i in 1:N.lambda) {

        # Récupération du résultat après application de la méthode SIR avec le lambda_i
        resSparseSIR <- SIR.threshold(Y, X, H = H, lambda = lambdas[i], thresholding = thresholding)

        # Stockage de l'estimation du vecteur beta dans la ligne i de la matrice
        mat.b.th[i,] <- resSparseSIR$beta

        # Stockage du nombre de 0 dans l'estimation du vecteur beta après seuillage pour
        # cette valeur de lambda
        vect.nb.zeros[i] <- resSparseSIR$nb.zeros

        # Stockage du cos carré entre b et b_seuillage pour cette valeur de lambda
        vect.cosca[i] <- resSparseSIR$cos.squared

        # Stocage des variables utiles pour cette valeur de lambda à l'indice i de la liste
        # des variables utiles
        list.relevant.variables[[i]] <- resSparseSIR$list.relevant.variables
    }


    # Création d'un vecteur qui contient p valeurs : pour chaque variable x_i
    # est associé l'indice du lambda à partir duquel la variable devient inutile. 
    # Cet indice correspond donc également au nombre de lambda pour lesquels la variable
    # est utile.
    indice.0 = colSums(mat.b.th / mat.b.th, na.rm = TRUE)

    # On recherche ensuite un point de rupture dans la liste indice.0 ordonnée. Ce point de rupture 
    # fixe le nombre de variables à ne pas sélectionner, soit le nombre de variable 
    # où l'estimation de beta donne zero pour un lambda donné.
    require(strucchange)
    fit_bp <- breakpoints(sort(indice.0, decreasing = FALSE) ~ 1, breaks = 1, h = 2 / p)

    # Si le nombre de coeficients de l'estimation de beta égaux à 0 et qui sont égaux au point
    # de rupture, donc au nombre de variable à ne pas sélectionner : 
    if (length(which(vect.nb.zeros == fit_bp$breakpoints)) > 0) {
        # l'indice optimal de lambda est l'indice minimum où le nombre de 0 apparaissant pour une valeur
        # de lambda est égal au nombre de valeurs inutiles au point de rupture
        indice.opt <- min(which(vect.nb.zeros == fit_bp$breakpoints))
    } else {
        # Sinon, on somme le nombre de fois où le nombre de 0 qui apparait dans l'estimation
        # de b pour chaque lambda est inférieur au nombre de valeurs inutiles au point de rupture
        # + 1
        indice.opt <- sum(vect.nb.zeros < fit_bp$breakpoints) + 1
    }

    # On récupère le lambda qui correspond à cet indice
    lambda.opt <- lambdas[indice.opt]

    # Affichage de la sélection du lambda
    if (graphic == TRUE) {

        if (length(which(vect.nb.zeros == fit_bp$breakpoints)) > 0) {
            par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
            # affichage des indices triés
            plot(sort(indice.0, decreasing = FALSE), xlab = "variable i", ylab = expression("Indices of" ~ lambda))
            # ligne verticale pour montrer la rupture
            abline(v = fit_bp$breakpoints + 0.5, col = 6, lwd = 3)
            title(paste("Choix de l'indice du lambda optimal parmi les", N.lambda, "testés"))
        }
        else {
            par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
            # affichage des indices triés
            plot(sort(indice.0, decreasing = FALSE), xlab = "variable p_i", ylab = "Indice du lambda à partir duquel la variable est inutile")
            # ligne verticale pour montrer la rupture
            abline(v = fit_bp$breakpoints + 1.5, col = 7, lwd = 3)
            title(paste("Choix de l'indice du lambda optimal parmi les ", N.lambda, " testés"))
        }

        # Affichage du pourcentage de variable utiles en fonction des lambdas
        plot(lambdas, 1 - vect.nb.zeros / p, ylim = c(0, 1.3), xlim = c(0, max(lambdas)), xlab = expression(lambda), col = 3, type = "l", ylab = "", lwd = 3)

        # Affichage de l'évolution du cos carré entre SIR et SParseSIR en fonction des valeurs de lambda
        points(lambdas, vect.cosca, type = "l", ylab = "", xlab = "lambda", lwd = 3, col = "black")
        title(paste("Evolution of cos^2 as a function of lambda (with", thresholding, "thresholding)"))

        # Ligne verticale au lambda optimal
        abline(v = lambda.opt, col = 6, lwd = 3)

        legend("topright", legend = c("% of relevant variables", expression("optimal" ~ lambda), expression(cos ^ 2 ~ (hat(b)[thresholding] ~ "," ~ hat(b)[SIR]))), col = c(3, 6, "black"), lty = c(1, 1, 1), lwd = c(3, 3, 3), cex = 1)
    }

    if (is.na(lambda.opt) == TRUE) {
        b.opt <- SIR(Y, X, H = 10)$beta

        list.relevant.var <- colnames(X)


    }
    else {
        # Cas où on trouve bien un lambda optimal parmi les N.lambda
        # Si le nombre de variables inutiles est inférieur au nombre de variable total-1
        if (fit_bp$breakpoints < (p - 1)) {

            # Le beta optimal est la ligne de la matrice d'estimations de beta qui 
            # correspond à lambda.opt
            b.opt <- mat.b.th[which(lambdas == lambda.opt),]
            # Conversion du beta en matrice à une ligne p colonnes
            b.opt <- matrix(b.opt, nrow = 1)
            # Renommage des colonnes
            colnames(b.opt) <- colnames(X)
            # Récupération des variables utiles qui sont les colonnes du b.opt
                       # où la valeur est différente de zero
            list.relevant.var <- colnames(b.opt)[-which(b.opt == 0)]

     
        }
        # Cas où la méthode Sparse SIR avec seuillage n'a pas permis de réaliser de la sélection de variable
        else {
            b.opt <- SIR(Y, X, H = 10)
        }
    }

    res = list(b.opt = b.opt, lambdas = lambdas, lambda.opt = lambda.opt, mat.b.th = mat.b.th,
               N.lambda = N.lambda, vect.nb.zeros = vect.nb.zeros,
               list.relevant.variables = list.relevant.var,n=n,p=p,H=H,M1=M1)

    class(res) = "SIR.threshold.opt"
    return(res)
}
