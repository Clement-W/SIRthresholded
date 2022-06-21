#' @export
#' @keywords internal
print.SIR = function(x, ...){
    if(!inherits(x,"SIR")) stop("Only use with \"SIR\" obects")
    
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    res = matrix(x$beta,ncol=1)
    rownames(res) = 1:x$p
    colnames(res) = "Estimated direction"
    
    cat("\n")
    prmatrix(signif(res,3))
    cat("\n")
}