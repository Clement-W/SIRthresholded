#' @export
#' @keywords internal
summary.SIR <- function(x, ...) {
    if (!inherits(x, "SIR")) stop("Only use with \"SIR\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")


    cat("\n===== Sliced Inverse Regression =====", "\n")

    cat("\n")
    cat(paste("Number of observations:", x$n), "\n")
    cat(paste("Number of variables in X:", x$p), "\n")
    cat(paste("Number of slices:", x$H), "\n")
    cat("\n")

    cat("Results of EDR directions estimation:\n")
    res <- matrix(x$b, ncol = 1)
    colnames(res) <- "Estimated direction"



    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")

}