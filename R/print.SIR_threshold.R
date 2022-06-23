#' @export
#' @keywords internal
print.SIR_threshold <- function(x, ...) {
    if (!inherits(x, "SIR_threshold"))
        stop("Only use with \"SIR_threshold\" obects")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    res <- matrix(x$b, ncol = 1)
    rownames(res) <- 1:x$p
    colnames(res) <- "Estimated direction"

    cat("\n")
    prmatrix(signif(res, 3))
    cat("\n")
}