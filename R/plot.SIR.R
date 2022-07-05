#'  Graphical output of SIR
#'
#' Display the 10 first eigen values and the estimated index versus Y of the SIR model
#' @param x A SIR object
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "eigvals" Plot the eigen values of the matrix of interest
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y
#'   \item "" Plot every graphs
#' }
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
#' res = SIR(Y, X, H = 10, graphic = FALSE)
#' 
#' # Eigen values
#' plot.SIR(res,choice="eigvals")
#'
#' # Estimated index versus Y
#' plot.SIR(res,choice="estim_ind")
#' @export
plot.SIR <- function(x, choice = "", ...) {

    if (!inherits(x, "SIR"))
        stop("Only use with \"SIR\" obects")

    if (!(choice %in% c("eigvals", "estim_ind", "")))
        stop("\"choice\" must be either \"eigvals\" or \"estim_ind\"", call. = FALSE)

    if (choice == "" || choice == "eigvals") {
        dev.new()


        eig_vals <- x$eig_val
        title <- "Eigenvalues"
        # Take only the 10 first eigen values
        if (length(eig_vals) >= 10) {
            eig_vals <- eig_vals[1:10]
            title <- "10 first eigenvalues"
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