#' @export
#' @keywords internal
summary.SIR <- function(object,...) {
    if (!inherits(object, "SIR")) stop("Only use with \"SIR\" obects")

    cat("\nCall:\n", deparse(object$call), "\n", sep = "")


    cat("\n===== Sliced Inverse Regression =====", "\n")

    cat("\n")
    cat(paste("Number of observations:", object$n), "\n")
    cat(paste("Number of variables in X:", object$p), "\n")
    cat(paste("Number of slices:", object$H), "\n")
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

}