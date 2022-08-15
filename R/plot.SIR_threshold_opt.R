#'  Graphical output of SIR_threshold_opt
#'
#' Display the 10 first eigen values,the estimated index versus Y of the SIR model,
#' the evolution of \eqn{cos^2} and variable selection according to \eqn{\lambda}, and the 
#' regularization path of \eqn{\hat{b}}.
#' @param x A SIR_threshold_opt object
#' @param choice the graph to plot: 
#' \itemize{
#'   \item "estim_ind" Plot the estimated index by the SIR model versus Y.
#'   \item "opt_lambda" Plot the choice of \eqn{\lambda_{opt}}.
#'   \item "cos2_selec" Plot the evolution of \eqn{cos^2} and variable selection according to 
#'   \eqn{\lambda}.
#'   \item "regul_path" Plot the regularization path of \eqn{\hat{b}}.
#'   \item "" Plot every graphs (default).
#' }
#' @param \ldots arguments to be passed to methods, such as graphical parameters (not used here).
#' @return No return value
#' @examples
#' # Generate Data
#' set.seed(10)
#' n <- 200
#' beta <- c(1,1,rep(0,8))
#' X <- mvtnorm::rmvnorm(n,sigma=diag(1,10))
#' eps <- rnorm(n)
#' Y <- (X%*%beta)**3+eps
#'
#' # Apply SIR with soft thresholding
#' res = SIR_threshold_opt(Y,X,H=10,n_lambda=100,thresholding="soft")
#' 
#' # Estimated index versus Y
#' plot(res,choice="estim_ind")
#' 
#' # Choice of optimal lambda
#' plot(res,choice="opt_lambda")
#' 
#' # Evolution of cos^2 and var selection according to lambda
#' plot(res,choice="cos2_selec")
#' 
#' # Regularization path
#' plot(res,choice="regul_path")
#' @export
#' @importFrom grDevices dev.new
#' @importFrom graphics abline axis box legend matplot mtext par points text title
plot.SIR_threshold_opt <- function(x, choice = "", ...) {

    if (!inherits(x, "SIR_threshold_opt"))
        stop("Only use with \"SIR_threshold_opt\" obects")

    if (!(choice %in% c("", "estim_ind", "opt_lambda", "cos2_selec", "regul_path")))
        stop("\"choice\" must be either \"estim_ind\",\"opt_lambda\",\"cos2_selec\" or
         \"regul_path\"", call. = FALSE)
    
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    if (choice == "" || choice == "estim_ind") {
        if (choice == "") {
            dev.new()
        }
        plot(x$index_pred, x$Y, xlab = "Estimated first index", ylab = "y", pch = 4)
        title("Reconstructed index")
    }

    if (choice == "" || choice == "opt_lambda") {
        if (choice == "") {
            dev.new()
        }
        par(las = 2)
        # Display sorted lambda's indices
        plot(sort(((x$indices_useless_var) / x$n_lambda) * 100, decreasing = FALSE),
                    xlab = "variables", ylab = expression("Proportion of" ~ lambda ~
                    "that select the variables (percent)"), pch = 16, axes = FALSE)
        axis(1, labels = names(sort(x$indices_useless_var)),
            at = 1:length(x$indices_useless_var))
        axis(2)
        box()

        # vertical line to show the break
        if (length(which(x$vect_nb_zeros == x$fit_bp$breakpoints)) > 0) {
            abline(v = x$fit_bp$breakpoints + 0.5, col = 6, lwd = 3)
        } else {
            abline(v = x$fit_bp$breakpoints + 1.5, col = 7, lwd = 3)
        }
        title(paste("Choosing the optimal lambda (with", x$n_lambda, "lambdas)"))
    }

    if (choice == "" || choice == "cos2_selec") {
        if (choice == "") {
            dev.new()
        }
        # Display the percentage of useful variables according to the lambdas
        plot(x$lambdas, 1 - x$vect_nb_zeros / x$p, ylim = c(0, 1.3),
            xlim = c(0, max(x$lambdas)), xlab = expression(lambda), col = 3, type = "l",
            ylab = "", lwd = 3)

        # Display the evolution of cos^2 between SIR and SIR_thresholded according to
        # lambda
        points(x$lambdas, x$vect_cos_squared, type = "l", ylab = "", xlab = "lambda",
            lwd = 3, col = "black")
        title(paste("Variable selection and cos^2 (with", x$thresholding, "thresholding)"))

        # Vertical line at optimal lambda
        abline(v = x$lambda_opt, col = 6, lwd = 3)

        legend("topright", legend = c("% of relevant variables",
            expression("optimal" ~ lambda), expression(cos ^ 2 ~ (hat(b)[thresholding] ~
            "," ~ hat(b)[SIR]))), col = c(3, 6, "black"), lty = c(1, 1, 1),
            lwd = c(3, 3, 3), cex = 1)
    }

    if (choice == "" || choice == "regul_path") {
        if (choice == "") {
            dev.new()
        }
        par(mar = c(5.1, 4.1, 6, 2.1))

        # Puts all the beta estimation in the same direction
        mat_b <- x$mat_b
        j0 <- which.max(abs(mat_b[1,]))
        for (i in 1:nrow(mat_b)) {
            if (mat_b[i, j0] != 0) {
                mat_b[i,] <- mat_b[i,] * sign(mat_b[i, j0])
            }
        }

        # Create interval and labels for the axis
        interval <- x$lambdas
        vect_nb_zeros_padded <- c(0, x$vect_nb_zeros)
        id_breaks <- which(vect_nb_zeros_padded[-1] -
            vect_nb_zeros_padded[-length(vect_nb_zeros_padded)] != 0)
        interval <- interval[id_breaks]
        lab <- x$p - x$vect_nb_zeros[id_breaks]

        matplot(x$lambdas, mat_b, type = "l", lty = 1, xlab = expression(lambda),
            xlim = c(-0.02, max(x$lambdas)), ylim = c(min(mat_b), 1.1 * max(mat_b)),
            ylab = "Value of the coefficients of b")
        title("Regularization path", line = 4)
        text(rep(0, ncol(mat_b)), mat_b[1,], colnames(x$b), pos = 2)
        axis(side = 3, labels = lab, at = interval)
        mtext("Number of selected variables", side = 3, line = 2.3)

        abline(v = x$lambda_opt, col = 6, lwd = 3)
        legend("topright", legend = expression("optimal" ~ lambda), col = 6,
            lty = 1, lwd = 3, cex = 1)

    }
}