#' @export
#' @keywords internal
print.SIR <- function(x,...) {
    if (!inherits(x, "SIR"))
        stop("Only use with \"SIR\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    res <- matrix(x$b, ncol = 1)
    if (!is.null(colnames(x$b))) {
        rownames(res) <- colnames(x$b)
    } else {
        rownames(res) <- paste("X", 1:x$p, sep = "")
    }
    colnames(res) <- "Estimated direction"

    cat("\nResults of EDR directions estimation:\n")
    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")
}