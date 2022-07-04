#' SIR optimally thresholded
#'
#' Apply a single-index SIR (Sliced Inverse Regression) on (X,Y) with H slices, 
#' thresholded by an optimal lambda parameter. The optimal lambda is found among n_lambda 
#' that threshold the interest matrix. Then, for each variable in X, we store how many
#' lambda selects this variable in a vector of size p. Thus, we find a breakpoint in 
# this sorted vector, which indicates the optimal lambda.
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices.
#' @param n_lambda The number of lambda to test. The N.lambda tested lambdas are 
#' between 0 and the maximum value of the interest matrix.
#' @param thresholding The thresholding method (choose between hard, soft)
#' @param graphic A boolean, set to TRUE to plot graphs 
#' @param output A boolean, set to TRUE to print informations
#' @param choice the graph to plot: 
#' \item{estim_ind}{Plot the estimated index by the SIR model versus Y}
#' \item{opt_lambda}{Plot the choice of the optimal lambda}
#' \item{cos2_selec}{Plot the evolution of cos^2 and variable selection according to 
#' lambda}
#' \item{regul_path}{Plot the regularization path of b}
#' \item{""}{Plot every graphs}
#' @return An object of class SIR_threshold_opt, with attributes:
#' \item{b}{This is the optimal estimated EDR direction, which is the principal 
#' eigenvector of the interest matrix.}
#' \item{lambdas}{A vector that contains the tested lambdas}
#' \item{lambda_opt}{The optimal lambda}
#' \item{mat_b}{A matrix of size p*n_lambda that contains an estimation of beta 
#' in the columns for each lambda}
#' \item{n_lambda}{The number of lambda tested}
#' \item{vect_nb_zeros}{The number of 0 in b for each lambda}
#' \item{list_relevant_variables}{A list that contains the variables selected by 
#' the model}
#' \item{fit_bp}{An object of class breakpoints from the strucchange package,
#' that contains informations about the breakpoint which allows to deduce the
#' optimal lambda.}
#' \item{indices_useless_var}{A vector that contains p values: each variable is 
#' associated with the number of lambda that selects this variable.}
#' \item{vect_nb_zeros}{The number of 0 in b for each lambda}
#' \item{vect_cos_squared}{A vector that contains for each lambda,
#' the cosine squared between vanilla SIR and SIR thresholded}
#' \item{Y}{The response vector.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' \item{H}{The chosen number of slices.}
#' \item{M1}{The interest matrix thresholded with the optimal lambda.}
#' \item{thresholding}{The thresholding method used}
#' \item{call}{Unevaluated call to the function.}
#' \item{X_reduced}{The X data restricted to the variables selected by the model.
#' It can be used to estimate a new SIR model on the relevant variables to improve
#' the estimation of b.}
#' \item{index_pred}{The index b'X estimated by SIR}
#' @examples 
#' # Generate Data
#' set.seed(10)
#' n <- 200
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with hard thresholding
#' SIR_threshold_opt(Y,X,H=10,n_lambda=300,thresholding="hard")
#' 
#' # Apply SIR with soft thresholding
#' SIR_threshold_opt(Y,X,H=10,n_lambda=300,thresholding="soft")
#' @export
#' @importFrom strucchange breakpoints
SIR_threshold_opt <- function(Y, X, H = 10, n_lambda = 100, thresholding = "hard",
    graphic = TRUE, output = TRUE, choice = "") {

    cl <- match.call()

    n <- nrow(X)
    p <- ncol(X)

    if (is.null(colnames(X))) {
        colnames(X) <- paste("X", 1:p, sep = "")
    }

    # Estimation of b and the interest matrix with the classic SIR method
    res_SIR <- SIR(Y, X, H = 10, graphic = FALSE)
    b_sir <- res_SIR$b
    M1 <- res_SIR$M1

    # Creation of a list of lambdas going from 0 to the maximum absolute value 
    # of M1, with a total of 100 values.
    lambdas <- seq(0, max(abs(M1)), length.out = n_lambda + 1)[-(n_lambda + 1)]

    # Initialization of a matrix of size n_lambda*p that will contain the estimation 
    # of b for each the lambda
    mat_b <- matrix(0, ncol = p, nrow = n_lambda)

    # Initialization of a vector of size n_lambda to receive the number of 
    # 0 found in b for each lambda 
    vect_nb_zeros <- rep(NA, n_lambda)

    # Initialization of a vector of size n_lambda to receive the cos^2
    # between b_sir and b_threshold_sir
    vect_cos_squared <- rep(NA, n_lambda)

    # Initialization of the list of useful variables
    list_relevant_variables <- list()

    # Application of the SIR method with thresholding, with the n_lambda values
    # of lambdas
    for (i in 1:n_lambda) {

        # Get the result of SIR thresholded with lambda_i
        res_SIR_th <- SIR_threshold(Y, X, H = H, lambda = lambdas[i],
            thresholding = thresholding, graphic = FALSE)

        # Store the corresponding b in row i of the matrix
        mat_b[i,] <- res_SIR_th$b

        # Store the number of 0 found in the corresponding b for this lambda value
        vect_nb_zeros[i] <- res_SIR_th$nb_zeros

        # Store the cos^2 between this b and b_sir
        vect_cos_squared[i] <- res_SIR_th$cos_squared

        # Store useful variables for this lambda value at index i of 
        # the list of useful variables
        list_relevant_variables[[i]] <- res_SIR_th$list_relevant_variables
    }


    # Creation of a vector which contains p values: each variable is associated 
    # with the index of the lambda from which the variable becomes useless. 
    # This index also corresponds to the number of lambda for which the variable
    # is useful.
    indices_useless_var <- colSums(mat_b / mat_b, na.rm = TRUE)
    names(indices_useless_var) <- colnames(X)

    # We then look for a breakpoint in the ordered index_useless_var vector. This
    # breakpoints corresponds to the number of useless variables. From this breakpoints,
    # the variables are more difficult to toggle to 0 with a threshold.
    fit_bp <- breakpoints(sort(indices_useless_var, decreasing = FALSE) ~ 1,
        breaks = 1, h = 2 / p)

    # If the number of useless variables associated to a lambda 
    # corresponds to the breakpoint :   
    if (length(which(vect_nb_zeros == fit_bp$breakpoints)) > 0) {

        # The index of the optimal lambda is the first lambda where the corresponding
        # number of zero in b is equal to the breakpoint, i.e the number of 
        # useless variables
        indice_opt <- min(which(vect_nb_zeros == fit_bp$breakpoints))
    } else {

        # Else, we sum the number of times where the number of 0 in b is less than
        # the number of useless variables. Then, the index of the optimal lambda
        # corresponds to this sum+1
        indice_opt <- sum(vect_nb_zeros < fit_bp$breakpoints) + 1
    }


    # Get the optimal lambda
    lambda_opt <- lambdas[indice_opt]

    # If the optimal lambda is NA, keep the result of vanilla SIR
    if (is.na(lambda_opt) == TRUE) {
        b <- b_sir
        list_relevant_var <- colnames(X)
        M1_th = M1

    }
    else {

        # If the optimal lambda is found and the number of useless variable is
        # less than the number of total variables-1
        if (fit_bp$breakpoints < (p - 1)) {

            # The optimal beta is the row of mat_b that corresponds to lambda.opt
            b <- mat_b[which(lambdas == lambda_opt),]
            # Convert b into a one-line matrix
            b <- matrix(b, nrow = 1)
            # Rename columns
            colnames(b) <- colnames(X)
            # Get the relevant variables (the columns of b where the value is not 0)
            list_relevant_var <- colnames(b)[-which(b == 0)]
            M1_th = SIR_threshold(Y, X, H = H, lambda = lambda_opt,
            thresholding = thresholding, graphic = FALSE)$M1
        }

        # If thresholded SIR could not make variable selection, keep the result of
        # vanilla SIR
        else {
            b <- b_sir
            list_relevant_var <- colnames(X)
            M1_th = M1
        }
    }

    # Create the X reduced variable by restricting X to the relevant variables.
    X_reduced <- X[, list_relevant_var, drop = FALSE]

    # Estimated index
    index_pred <- X %*% t(b)


    res <- list(b = b, lambdas = lambdas, lambda_opt = lambda_opt,
        mat_b = mat_b, n_lambda = n_lambda, vect_nb_zeros = vect_nb_zeros,
        fit_bp = fit_bp, indices_useless_var = indices_useless_var,
        vect_cos_squared = vect_cos_squared, list_relevant_variables = list_relevant_var,
        n = n, p = p, H = H, M1 = M1_th, thresholding = thresholding, call = cl,
        X_reduced = X_reduced, index_pred = index_pred, Y = Y)

    class(res) <- "SIR_threshold_opt"

    if (graphic == TRUE) {
        plot.SIR_threshold_opt.R(res, choice = choice)
    }

    return(res)
}
