#' @export
#' @keywords internal
print.SIR.threshold.bootstrap = function(x, ...){
    if(!inherits(x,"SIR.threshold.bootstrap")) stop("Only use with \"SIR.threshold.bootstrap\" obects")
    
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    
    res = matrix(x$b.opt,ncol=1)
    rownames(res) = 1:x$p
    colnames(res) = "Estimated direction"
    
    cat("\n")
    prmatrix(signif(res,3))
    cat("\n")
}