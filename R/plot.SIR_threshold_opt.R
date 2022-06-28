#' @export
#' @keywords internal
plot.SIR_threshold_opt.R <- function(x, ...) {

    if (!inherits(x, "SIR_threshold_opt"))
        stop("Only use with \"SIR_threshold_opt\" obects")


    dev.new()
    plot(x$index_pred, x$Y, xlab = "Estimated index", ylab = "y", pch = 4)
    title("Reconstructed index")


    if (length(which(x$vect_nb_zeros == x$fit_bp$breakpoints)) > 0) {
        dev.new()
        # affichage des indices triés
        plot(sort(((x$indices_useless_var) / x$n_lambda) * 100, decreasing = FALSE), xlab = "variable i",
                ylab = expression("Proportion of" ~ lambda ~
                "that select the variable i (percent)"), pch = 16)

        # ligne verticale pour montrer la rupture
        abline(v = x$fit_bp$breakpoints + 0.5, col = 6, lwd = 3)
        title(paste("Choosing the optimal lambda (with", x$n_lambda, "lambdas)"))


    } else {
        dev.new()
        # affichage des indices triés
        plot(sort(((x$indices_useless_var) / x$n_lambda) * 100, decreasing = FALSE), xlab = "variable i",
             ylab = expression("Proportion of" ~ lambda ~
             "that select the variable i (percent)"))
        # ligne verticale pour montrer la rupture
        abline(v = x$fit_bp$breakpoints + 1.5, col = 7, lwd = 3)
        title(paste("Choosing the optimal lambda (with", x$n_lambda, "lambdas)"))
    }

    dev.new()
    # Affichage du pourcentage de variable utiles en fonction des lambdas
    plot(x$lambdas, 1 - x$vect_nb_zeros / x$p, ylim = c(0, 1.3),
        xlim = c(0, max(x$lambdas)), xlab = expression(lambda), col = 3, type = "l",
        ylab = "", lwd = 3)

    # Affichage de l'évolution du cos carré entre SIR et SParseSIR en fonction de lambda
    points(x$lambdas, x$vect_cosca, type = "l", ylab = "", xlab = "lambda", lwd = 3,
        col = "black")
    title(paste("Variable selection and cos² (with", x$thresholding, "thresholding)"))

    # Ligne verticale au lambda optimal
    abline(v = x$lambda_opt, col = 6, lwd = 3)

    legend("topright", legend = c("% of relevant variables",
        expression("optimal" ~ lambda), expression(cos ^ 2 ~ (hat(b)[thresholding] ~
        "," ~ hat(b)[SIR]))), col = c(3, 6, "black"), lty = c(1, 1, 1),
        lwd = c(3, 3, 3), cex = 1)


    dev.new()
    par(mar = c(5.1, 4.1, 6, 2.1))
    mat_b <- x$mat_b
    j0 <- which.max(abs(mat_b[1,]))
    for (i in 1:nrow(mat_b)) {
        mat_b[i,] = mat_b[i,] * sign(mat_b[i, j0])
    }

    interval = x$lambdas
    vect_nb_zeros_padded = c(0, x$vect_nb_zeros)
    id_ruptures = which(vect_nb_zeros_padded[-1] - vect_nb_zeros_padded[-length(vect_nb_zeros_padded)] != 0)
    interval = interval[id_ruptures]
    lab = x$p - x$vect_nb_zeros[id_ruptures]

    matplot(x$lambdas, mat_b, type = "l", lty = 1, xlab = expression(lambda), xlim = c(-0.02, max(x$lambdas)), ylim = c(0, 1.1))
    title("Regularization path", line = 4)
    text(rep(0, ncol(mat_b)), mat_b[1,], colnames(x$b), pos = 2)
    axis(side = 3, labels = lab, at = interval)
    mtext("Number of variables selected", side = 3, line = 2.3)

    abline(v = x$lambda_opt, col = 6, lwd = 3)
    legend("topright", legend = expression("optimal" ~ lambda), col = 6, lty = 1, lwd = 3, cex = 1)
}