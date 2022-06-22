#' @export
#' @keywords internal
plot.SIR.threshold.opt.R <- function(x, ...) {

    if (!inherits(x, "SIR.threshold.opt"))
        stop("Only use with \"SIR.threshold.opt\" obects")


    if (length(which(x$vect.nb.zeros == x$fit_bp$breakpoints)) > 0) {
        par(mfrow = c(1, 2))
        # affichage des indices triés
        plot(sort(x$indice.0, decreasing = FALSE), xlab = "variable i",
                ylab = expression("Indices of" ~ lambda ~
                "from which the variable is not selected"))

        # ligne verticale pour montrer la rupture
        abline(v = x$fit_bp$breakpoints + 0.5, col = 6, lwd = 3)
        title(paste("Choosing the optimal lambda (with", x$N.lambda, "lambdas)"))

    } else {
        par(mfrow = c(1, 2))
        # affichage des indices triés
        plot(sort(x$indice.0, decreasing = FALSE), xlab = "variable p_i",
            ylab = expression("Indices of" ~ lambda ~
            "from which the variable i is not selected"))
        # ligne verticale pour montrer la rupture
        abline(v = x$fit_bp$breakpoints + 1.5, col = 7, lwd = 3)
        title(paste("Choosing the optimal lambda (with", x$N.lambda, "lambdas)"))
    }

    # Affichage du pourcentage de variable utiles en fonction des lambdas
    plot(x$lambdas, 1 - x$vect.nb.zeros / x$p, ylim = c(0, 1.3),
        xlim = c(0, max(x$lambdas)), xlab = expression(lambda), col = 3, type = "l",
        ylab = "", lwd = 3)

    # Affichage de l'évolution du cos carré entre SIR et SParseSIR en fonction de lambda
    points(x$lambdas, x$vect.cosca, type = "l", ylab = "", xlab = "lambda", lwd = 3,
        col = "black")
    title(paste("Variable selection and cos² (with", x$thresholding, "thresholding)"))

    # Ligne verticale au lambda optimal
    abline(v = x$lambda.opt, col = 6, lwd = 3)

    legend("topright", legend = c("% of relevant variables",
        expression("optimal" ~ lambda), expression(cos ^ 2 ~ (hat(b)[thresholding] ~
        "," ~ hat(b)[SIR]))), col = c(3, 6, "black"), lty = c(1, 1, 1),
        lwd = c(3, 3, 3), cex = 1)
}