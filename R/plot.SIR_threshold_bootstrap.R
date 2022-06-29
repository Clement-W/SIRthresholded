#' @export
#' @keywords internal
plot.SIR_threshold_bootstrap <- function(x, ...) {

    if (!inherits(x, "SIR_threshold_bootstrap"))
        stop("Only use with \"SIR_threshold_bootstrap\" obects")

    # Histogramme du nombre de variables sélectionnées par le modèle
    dev.new()
    barplot((table(x$nb_var_selec) / x$n_replications) * 100, ylab = "percent",
            xlab = "Number of selected variables")
    title("Sizes of the reduced models")

    dev.new()
    # Histogramme du nombre de fois ou chaque variable a été sélectionnée
    barplot((x$effectif_var / x$n_replications) * 100, names.arg =
            colnames(x$b), ylab = "percent", xlab = "variable name")
    title("Selected variables in the reduced models")
    
    dev.new()
    mat_b <- x$mat_b
    j0 <- which.max(abs(mat_b[1,]))
    for (i in 1:nrow(mat_b)) {
        mat_b[i,] = mat_b[i,] * sign(mat_b[i, j0])
    }
    boxplot(mat_b,xlab="coefficients of b",ylab="",main=paste("Value of b over",x$n_replications,"replications"),names=colnames(x$b))
    
    dev.new()
    boxplot(x$lambdas_opt_boot,xlab=expression(lambda),ylab="",main=paste("Value of optimal" ,expression(lambda), "over the",x$n_replications,"replications"))

}