#' @export
#' @keywords internal
print.SIR.threshold <- function(x, ...) {
    if (!inherits(x, "SIR.threshold"))
        stop("Only use with \"SIR.threshold\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    res <- matrix(x$b, ncol = 1)
    rownames(res) <- 1:x$p
    colnames(res) <- "Estimated direction"

    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")
}