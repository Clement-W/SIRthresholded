#' @export
#' @keywords internal
summary.SIR.threshold.bootstrap = function(x, ...){
    if(!inherits(x,"SIR.threshold.bootstrap")) stop("Only use with \"SIR.threshold.bootstrap\" obects")
    
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    
    cat("\n===== Optimally Thresholded Sliced Inverse Regression on bootstrapped replications =====","\n")
    
    cat("\n")
    cat(paste("Number of observations:",x$n),"\n")
    cat(paste("Number of variables in X:",x$p),"\n")
    cat(paste("Number of slices:",x$H),"\n")
    cat(paste("Number of bootstraped replications:",x$Nb.replications,"\n"))
    cat(paste("Optimal thresholding parameter lambda :",x$lambda.opt),"\n")
    cat(paste("Thresholding method:",x$thresholding))
    cat("\n\n")
    
    cat(paste("Number of selected variables = ", x$Nb.var.selec.opt, " over the ", x$p, " available variables"), fill = TRUE)
    cat(paste("List of relevant variables:", paste(x$list.relevant.variables, collapse = ",")),"\n")
    cat("\n")
    
    cat("Results of EDR directions estimation:\n")
    res = matrix(x$b.opt,ncol=1)
    rownames(res) = 1:x$p
    colnames(res) = "Estimated direction"
    
    cat("\n")
    prmatrix(signif(res,3))
    cat("\n")
}