#'  Graphical output of SIR_bootstrap
#'
#' Display the 10 first eigen values and the estimated index versus Y of the 
#' SIRbootstrap model.
#' @param x A SIR_bootstrap object
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "eigvals" Plot the eigen values of the matrix of interest.
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y.
#'   \item "" Plot every graphs (default).
#' }
#' @param \ldots arguments to be passed to methods, such as graphical parameters (not used here).
#' @return No return value
#' @examples
#' # Generate Data
#' set.seed(10)
#' n <- 500
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply bootstrap SIR
#' res = SIR_bootstrap(Y, X, H = 10, B = 10)
#' 
#' # Eigen values
#' plot(res,choice="eigvals")
#'
#' # Estimated index versus Y
#' plot(res,choice="estim_ind")
#' @export
#' @importFrom grDevices dev.new
#' @importFrom graphics axis barplot title
plot.SIR_bootstrap <- function(x, choice = "", ...) {

    if (!inherits(x, "SIR_bootstrap"))
        stop("Only use with \"SIR_bootstrap\" obects")

    if (!(choice %in% c("eigvals", "estim_ind", "")))
        stop("\"choice\" must be either \"eigvals\" or \"estim_ind\"", call. = FALSE)

    if (choice == "" || choice == "eigvals") {
        dev.new()

        eig_vals <- x$eig_val
        title <- "Eigenvalues"
        # Take only the 10 first eigen values
        if (length(eig_vals) >= 10) {
            eig_vals <- eig_vals[1:10]
            title <- "10 largest eigenvalues"
        }

        a <- barplot(eig_vals, ylab = "eigenvalues", xlab = "dimensions", main = title)
        axis(1, at = a, labels = seq(1, length(eig_vals)))
    }


    if (choice == "" || choice == "estim_ind") {
        dev.new()
        plot(x$index_pred, x$Y, xlab = "Estimated first index", ylab = "y", pch = 4)
        title("Reconstructed index")
    }

}