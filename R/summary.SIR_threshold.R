#' @export
#' @keywords internal
summary.SIR_threshold = function(x, ...) {
    if (!inherits(x, "SIR_threshold"))
        stop("Only use with \"SIR_threshold\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    cat("\n===== Thresholded Sliced Inverse Regression =====", "\n")

    cat("\n")
    cat(paste("Number of observations:", x$n), "\n")
    cat(paste("Number of variables in X:", x$p), "\n")
    cat(paste("Number of slices:", x$H), "\n")
    cat(paste("Thresholding parameter lambda:", signif(x$lambda, 3)), "\n")
    cat(paste("Thresholding method:", x$thresholding))
    cat("\n")

    cat("\n")
    cat(paste("Number of selected variables = ", length(x$list_relevant_variables),
        " over the ", x$p, " available variables"), fill = TRUE)
    cat(paste("List of relevant variables:", paste(x$list_relevant_variables,
        collapse = ",")), "\n")
    cat("\n")

    cat("Results of EDR directions estimation:\n")
    res <- matrix(x$b, ncol = 1)
    if (!is.null(colnames(x$b))) {
        rownames(res) <- colnames(x$b)
    } else {
        rownames(res) <- paste("X", 1:x$p, sep = "")
    }
    colnames(res) <- "Estimated direction"

    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")

}