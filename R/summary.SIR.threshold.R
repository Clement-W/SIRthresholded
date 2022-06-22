#' @export
#' @keywords internal
summary.SIR.threshold = function(x, ...) {
    if (!inherits(x, "SIR.threshold"))
        stop("Only use with \"SIR.threshold\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    cat("\n===== Thresholded Sliced Inverse Regression =====", "\n")

    cat("\n")
    cat(paste("Number of observations:", x$n), "\n")
    cat(paste("Number of variables in X:", x$p), "\n")
    cat(paste("Number of slices:", x$H), "\n")
    cat(paste("Thresholding parameter lambda:", x$lambda), "\n")
    cat(paste("Thresholding method:", x$thresholding))
    cat("\n")

    cat("\n")
    cat(paste("Number of selected variables = ", length(x$list.relevant.variables),
        " over the ", x$p, " available variables"), fill = TRUE)
    cat(paste("List of relevant variables:", paste(x$list.relevant.variables,
        collapse = ",")), "\n")
    cat("\n")

    cat("Results of EDR directions estimation:\n")
    res <- matrix(x$beta, ncol = 1)
    rownames(res) <- 1:x$p
    colnames(res) <- "Estimated direction"

    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")

}