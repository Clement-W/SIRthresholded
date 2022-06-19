#' @export
#' @keywords internal
print.SIR = function(x, ...){
    if(!inherits(x,"SIR")) stop("Only use with \"SIR\" obects")
    
    cat("\n===== Perform Sliced Inverse Regression =====","\n")
    
    cat("\n")
    cat(paste("Number of observations:",x$n),"\n")
    cat(paste("Number of variables in X:",x$p),"\n")
    cat(paste("Number of slices:",x$H),"\n")
    cat("\n")
    
    cat("Results of EDR directions estimation:\n")
    res = matrix(x$beta,ncol=1)
    rownames(res) = 1:x$p
    colnames(res) = "Estimated direction"
    
    cat("\n")
    prmatrix(signif(res,3))
    cat("\n")
    
}
# Ça fonctionne bien mais le fait que la fonction soit pas load dans le global environement ça ne marche pas
