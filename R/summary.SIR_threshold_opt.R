#' @export
#' @keywords internal
summary.SIR_threshold_opt <- function(x, ...) {
    if (!inherits(x, "SIR_threshold_opt"))
        stop("Only use with \"SIR_threshold_opt\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    cat("\n===== Optimally Thresholded Sliced Inverse Regression =====", "\n")

    cat("\n")
    cat(paste("Number of observations:", x$n), "\n")
    cat(paste("Number of variables in X:", x$p), "\n")
    cat(paste("Number of slices:", x$H), "\n")
    cat(paste("Testing:", x$n_lambda, "lambda between 0 and", max(abs(x$M1))), "\n")
    cat(paste("Optimal thresholding parameter lambda :", x$lambda_opt), "\n")
    cat(paste("Thresholding method:", x$thresholding, "\n"))
    cat("\n")

    cat(paste("Number of selected variables = ", length(x$list_relevant_variables),
        " over the ", x$p, " available variables"), fill = TRUE)
    cat(paste("List of relevant variables:", paste(x$list_relevant_variables,
        collapse = ",")), "\n")
    cat("\n")

    cat("Results of EDR directions estimation:\n")
    res <- matrix(x$b, ncol = 1)
    rownames(res) <- 1:x$p
    colnames(res) <- "Estimated direction"

    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")
}