#' @export
#' @keywords internal
summary.SIR.threshold.opt = function(x, ...){
    if(!inherits(x,"SIR.threshold.opt")) stop("Only use with \"SIR.threshold.opt\" obects")
    
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    
    cat("\n===== Optimally Thresholded Sliced Inverse Regression =====","\n")
    
    cat("\n")
    cat(paste("Number of observations:",x$n),"\n")
    cat(paste("Number of variables in X:",x$p),"\n")
    cat(paste("Number of slices:",x$H),"\n")
    cat(paste("Testing:",x$N.lambda,"lambda between 0 and",max(abs(x$M1))),"\n")
    cat(paste("Optimal thresholding parameter lambda :",x$lambda.opt),"\n")
    cat(paste("Thresholding method:",x$thresholding,"\n"))
    cat("\n")
    
    cat(paste("Number of selected variables = ", length(x$list.relevant.variables), " over the ", x$p, " available variables"), fill = TRUE)
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