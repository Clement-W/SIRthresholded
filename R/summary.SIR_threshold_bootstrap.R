#' @export
#' @keywords internal
summary.SIR_threshold_bootstrap <- function(object, ...) {
    if (!inherits(object, "SIR_threshold_bootstrap"))
        stop("Only use with \"SIR_threshold_bootstrap\" obects")

    cat("\nCall:\n", deparse(object$call), "\n", sep = "")

    cat(paste("\n===== Optimally Thresholded Sliced Inverse Regression on",
        "bootstrapped replications =====", "\n"))

    cat("\n")
    cat(paste("Number of observations:", object$n), "\n")
    cat(paste("Number of variables in X:", object$p), "\n")
    cat(paste("Number of slices:", object$H), "\n")
    cat(paste("Number of bootstraped replications:", object$n_replications, "\n"))
    cat(paste("Optimal thresholding parameter lambda :", signif(object$lambda_opt, 3)), "\n")
    cat(paste("Thresholding method:", object$thresholding))
    cat("\n\n")

    cat(paste("Number of selected variables = ", object$nb_var_selec_opt, " over the ",
        object$p, " available variables"), fill = TRUE)
    cat(paste("List of relevant variables:", paste(object$list_relevant_variables,
        collapse = ",")), "\n")
    cat("\n")

    cat("Results of EDR directions estimation:\n")
    res <- matrix(object$b, ncol = 1)
    if (!is.null(colnames(object$b))) {
        rownames(res) <- colnames(object$b)
    } else {
        rownames(res) <- paste("X", 1:object$p, sep = "")
    }
    colnames(res) <- "Estimated direction"

    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")

    cat("Estimate a new SIR model on the relevant variables with :\n")
    varname = deparse(substitute(object))
    cat(paste("\t SIR(Y=", varname, "$Y, X=", varname, "$X_reduced, H=", varname, "$H)\n", sep = ""))
}