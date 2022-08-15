#' SIR optimally thresholded on bootstraped replications
#'
#' Apply a single-index optimally soft/hard thresholded \eqn{SIR} with \eqn{H} slices on 
#' `n_replications` bootstraped replications of \eqn{(X,Y)}. The optimal number of 
#' selected variables is the number of selected variables that came back most often 
#' among the replications performed. From this, we can get the corresponding \eqn{\hat{b}} 
#' and \eqn{\lambda_{opt}} that produce the same number of selected variables in the result of 
#' `SIR_threshold_opt`. 
#' @param X A matrix representing the quantitative explanatory variables (bind by column).
#' @param Y A numeric vector representing the dependent variable (a response vector).
#' @param H The chosen number of slices (default is 10).
#' @param n_lambda The number of lambda to test. The n_lambda tested lambdas are
#' uniformally distributed between 0 and the maximum value of the interest matrix (default is 100).
#' @param thresholding The thresholding method to choose between hard and soft (default is hard).
#' @param n_replications The number of bootstraped replications of (X,Y) done to
#' estimate the model (default is 50).
#' @param k Multiplication factor of the bootstrapped sample size
#' (default is 1 = keep the same size as original data).
#' @param graph A boolean, set to TRUE to plot graphs (default is TRUE).
#' @param output A boolean, set to TRUE to print information (default is TRUE).
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y.
#'   \item "size" Plot the size of the models across the replications.
#'   \item "selec_var" Plot the occurrence of the selected variables across the replications.
#'   \item "coefs_b" Plot the value of b across the replications.
#'   \item "lambdas_replic" Plot the optimal lambdas across the replications.
#'   \item "" Plot every graphs (default).
#' }
#' @return An object of class SIR_threshold_bootstrap, with attributes:
#' \item{b}{This is the optimal estimated EDR direction, which is the principal 
#' eigenvector of the interest matrix.}
#' \item{lambda_opt}{The optimal lambda.}
#' \item{vec_nb_var_selec}{Vector that contains the number of selected variables 
#' for each replications.}
#' \item{occurrences_var}{Vector that contains at index i the number of times the 
#' i_th variable has been selected in a replication.}
#' \item{call}{Unevaluated call to the function.}
#' \item{nb_var_selec_opt}{Optimal number of selected variables which is the number 
#' of selected variables that came back most often among the replications performed.}
#' \item{list_relevant_variables}{A list that contains the variables selected by 
#' the model.}
#' \item{n}{Sample size.}
#' \item{p}{The number of variables in X.}
#' \item{H}{The chosen number of slices.}
#' \item{n_replications}{The number of bootstraped replications of (X,Y) done to
#' estimate the model.}
#' \item{thresholding}{The thresholding method used.}
#' \item{X_reduced}{The X data restricted to the variables selected by the model.
#' It can be used to estimate a new SIR model on the relevant variables to improve
#' the estimation of b.}
#' \item{mat_b}{Contains the estimation b at each bootstraped replications.}
#' \item{lambdas_opt_boot}{Contains the optimal lambda found by SIR_threshold_opt at 
#' each replication.}
#' \item{index_pred}{The index Xb' estimated by SIR.}
#' \item{Y}{The response vector.}
#' \item{M1}{The interest matrix thresholded with the optimal lambda.}
#' @examples
#'
#' # Generate Data
#' set.seed(8)
#' n <-  170
#' beta <- c(1,1,1,1,1,rep(0,15))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,20))
#' eps <- rnorm(n,sd=8)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with hard thresholding
#' \donttest{SIR_threshold_bootstrap(Y,X,H=10,n_lambda=300,thresholding="hard", n_replications=30,k=2)}
#' 
#' @export
SIR_threshold_bootstrap <- function(Y, X, H = 10, thresholding = "hard",
    n_replications = 50, graph = TRUE, output = TRUE, n_lambda = 100, k = 2,
    choice = "") {

    cl <- match.call()

    # Ensure that X and Y are matrices
    X = ensure_matrix(X)
    Y = ensure_matrix(Y)

    # SIR optimally thresholded
    res_SIR_th <- SIR_threshold_opt(Y, X, H = H, thresholding = thresholding,
        graph = FALSE, output = FALSE, n_lambda = n_lambda)

    p <- ncol(X)
    n <- nrow(X)

    if (is.null(colnames(X))) {
        colnames(X) <- paste("X", 1:p, sep = "")
    }

    # Vector used to store the number of selected variables for each replications
    vec_nb_var_selec <- rep(NA, n_replications)

    # Initialize the list of relevant variables to every variables at first
    list_relevant_variables <- colnames(X)

    # List used to store the selected variables for each replications
    liste_var_selec <- list()

    # Vector that contains the number of times the variable at index j has been selected
    # in a replication
    occurrences_var <- rep(0, p)

    # To store the estimation b at each bootstraped replications
    mat_b <- matrix(0, ncol = p, nrow = n_replications)

    # To store the optimal lambda at each replication 
    lambdas_opt_boot <- rep(0, n_replications)


    # For each replication
    for (replic in 1:n_replications) {
        # bootstrap 
        indice_bootstrap <- sample(1:n, size = k * n, replace = TRUE)
        # Get bootstraped X and Y
        X_boot <- X[indice_bootstrap,]
        Y_boot <- Y[indice_bootstrap]

        # SIR optimally thresholded on the bootstraped sample
        res_boot <- SIR_threshold_opt(Y_boot, X_boot, H = H, thresholding = thresholding,
            graph = FALSE, output = FALSE, n_lambda = n_lambda)

        # Store the number of selected variables 
        vec_nb_var_selec[replic] <- length(res_boot$list_relevant_variables)

        # Store the selected variables
        liste_var_selec[[replic]] <- res_boot$list_relevant_variables

        # Store the vector b 
        mat_b[replic,] <- res_boot$b

        # Store the optimal lambda
        lambdas_opt_boot[replic] <- res_boot$lambda_opt

        if (output && replic %% 5 == 0) {
            print(paste("Replication", replic, "/", n_replications))
        }
    }

    # For each replication
    for (i in 1:n_replications) {
        # Pour eeach variable
        for (j in 1:p) {
            # If the variable at index j was selected on the i_th replication
            if ((sum(liste_var_selec[[i]] == colnames(X)[j])) == 1) {
                # Increment the occurrence counter of the variable j
                occurrences_var[j] <- occurrences_var[j] + 1
            }
        }
    }


    # Optimal number of selected variables = the number of selected variables 
    # that came back most often among the replications performed
    nb_var_selec_opt <-
        as.numeric(names(which(table(vec_nb_var_selec) == max(table(vec_nb_var_selec)))))

    # Number of optimal zero
    nb_zeros_opt <- p - nb_var_selec_opt

    # estimation of b by taking the beta estimated on the whole sample by the
    # method, at the lambda from which the optimal number of zeros appears
    if (sum(res_SIR_th$vect_nb_zeros == nb_zeros_opt) != 0) {
        b <- res_SIR_th$mat_b[min(which(res_SIR_th$vect_nb_zeros == nb_zeros_opt)),]
    }
    # If the nb_zero_opt did not appear for any lambda in res_SIR_th$vect_nb_zeros,
    # decrease nb_zero_opt to select less variable, until it exists in res_SIR_th$vect_nb_zeros
    else {
        while (sum(res_SIR_th$vect_nb_zeros == nb_zeros_opt) == 0) {
            nb_zeros_opt = nb_zeros_opt - 1
        }
        b <- res_SIR_th$mat_b[min(which(res_SIR_th$vect_nb_zeros == nb_zeros_opt)),]

    }
    # Convert it into a one-line matrix
    b <- matrix(b, nrow = 1)
    # Rename columns
    colnames(b) <- colnames(X)


    # If variables have been removed
    if (sum(nb_zeros_opt > 0) > 0) {
        # Update the list of relevant variables
        list_relevant_variables <- colnames(X)[-which(b == 0)]

        # If the number of zero in b is 0
        if (length(which(b == 0)) == 0) {
            # The list of relevant variables contains every variables
            list_relevant_variables <- colnames(X)
        }
    }

    # Create the X reduced variable by restricting X to the relevant variables.
    X_reduced <- X[, list_relevant_variables, drop = FALSE]

    # Estimated index
    index_pred <- X %*% t(b)

    # Get the optimal lambda that corresponds to the optimal number of zero found earlier
    lambda_optim <-
        res_SIR_th$lambdas[min(which(res_SIR_th$vect_nb_zeros == nb_zeros_opt))]

    # Get the corresponding thresholded interest matrix
    M1_th = SIR_threshold(Y, X, H = H, lambda = lambda_optim,
            thresholding = thresholding, graph = FALSE)$M1

    res <- list(b = b, lambda_opt = lambda_optim,
        vec_nb_var_selec = vec_nb_var_selec, occurrences_var = occurrences_var, call = cl,
        nb_var_selec_opt = nb_var_selec_opt, list_relevant_variables =
        list_relevant_variables, n = n, p = p, H = H, n_replications =
        n_replications, thresholding = thresholding, X_reduced = X_reduced, mat_b = mat_b,
        lambdas_opt_boot = lambdas_opt_boot, index_pred = index_pred, Y = Y, M1 = M1_th)

    class(res) <- "SIR_threshold_bootstrap"

    if (graph == TRUE) {
        plot.SIR_threshold_bootstrap(res, choice = choice)
    }

    return(res)
}
