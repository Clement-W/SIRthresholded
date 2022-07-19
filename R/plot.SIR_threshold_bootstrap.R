#'  Graphical output of SIR_threshold_bootstrap
#'
#' Display the estimated index versus Y of the SIR model, the size of the models,
#' the occurrence of variable selection, the distribution of the coefficients of  
#' and \eqn{\hat{b}} and the distribution of \eqn{\lambda_{opt}} found across the replications.
#' @param x A SIR_threshold_bootstrap object
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y
#'   \item "size" Plot the size of the models across the replications
#'   \item "selec_var" Plot the occurrence of the selected variables across the replications
#'   \item "coefs_b" Plot the value of \eqn{\hat{b}} across the replications
#'   \item "lambdas_replic" Plot the distribution of \eqn{\lambda_{opt}} across the replications
#'   \item "" Plot every graphs
#' }
#' @param \ldots arguments to be passed to methods, such as graphical parameters (not used here).
#' @examples
#' ## Generate Data
#' # set.seed(10)
#' # n <- 200
#' # beta <- c(1,1,rep(0,8))
#' # X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' # eps <- rnorm(n)
#' # Y <- (X%*%beta)**3+eps
#'
#' # res = SIR_threshold_bootstrap(Y,X,H=10,n_lambda=300,thresholding="hard", n_replications=30,k=2)
#'
#' ## Estimated index versus Y
#' # plot(res,choice="estim_ind")
#' 
#' ## Model size
#' # plot(res,choice="size")
#' 
#' ## Selected variables
#' # plot(res,choice="selec_var")
#' 
#' ## Coefficients of b
#' # plot(res,choice="coefs_b")
#'
#' ## Optimal lambdas
#' # plot(res,choice="lambdas_replic")
#' @export
#' @importFrom grDevices dev.new
#' @importFrom graphics abline barplot box boxplot legend title
plot.SIR_threshold_bootstrap <- function(x, choice = "",...) {

    if (!inherits(x, "SIR_threshold_bootstrap"))
        stop("Only use with \"SIR_threshold_bootstrap\" obects")

    if (!(choice %in% c("", "estim_ind", "size", "selec_var", "coefs_b",
        "lambdas_replic")))
        stop("\"choice\" must be either \"estim_ind\",\"size\",\"selec_var\",
            \"coefs_b\" or \"lambdas_replic\"", call. = FALSE)

    if (choice == "" || choice == "estim_ind") {
        if (choice == "") {
            dev.new()
        }
        plot(x$index_pred, x$Y, xlab = "Estimated first index", ylab = "y", pch = 4)
        title("Reconstructed index")
    }

    if (choice == "" || choice == "size") {
        # Barplot of the number of selected variable by the model 
        if (choice == "") {
            dev.new()
        }
        barplot((table(x$vec_nb_var_selec) / x$n_replications) * 100, ylab = "percent",
                xlab = "Number of selected variables")
        title(paste("Number of selected variables over the", x$n_replications, "bootstrap replications"))
    }

    if (choice == "" || choice == "selec_var") {
        if (choice == "") {
            dev.new()
        }
        # Barplot of the number of time where each variable has been selected
        barplot((x$occurrences_var / x$n_replications) * 100, names.arg =
                colnames(x$b), ylab = "percent", xlab = "variable name")
        title(paste("Selected variables over the", x$n_replications, "bootstrap replications"))
    }

    if (choice == "" || choice == "coefs_b") {
        # Boxplot showing the distribution of the coefficients of b
        if (choice == "") {
            dev.new()
        }
        mat_b <- x$mat_b
        j0 <- which.max(abs(mat_b[1,]))
        for (i in 1:nrow(mat_b)) {
            if(mat_b[i, j0] != 0){
                mat_b[i,] <- mat_b[i,] * sign(mat_b[i, j0])
            }
        }
        boxplot(mat_b, xlab = "coefficients of b", ylab = "", main =
            paste("Value of b over the", x$n_replications, "bootstrap replications"),
            names = colnames(x$b))
        points(matrix(x$b * sign(x$b[j0]), ncol = 1), pch = 19, col = 6, lwd = 3)
        legend("topright", legend = expression(hat(beta) ~ " associated to optimal" ~ lambda), col = 6,
               pch = 19)
    }

    if (choice == "" || choice == "lambdas_replic") {
        # Boxplot showing the distribution of the optimal lambdas found
        if (choice == "") {
            dev.new()
        }
        boxplot(x$lambdas_opt_boot, xlab = expression(lambda), ylab = "",
            main = paste("Value of optimal", expression(lambda), "over the",
            x$n_replications, "bootstrap replications"))
        abline(h = x$lambda_opt, col = 6, lwd = 3)
        legend("topright", legend = expression("optimal" ~ lambda), col = 6,
               lty = 1, lwd = 3, cex = 1)
    }
}