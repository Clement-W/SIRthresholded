#' @export
#' @keywords internal
plot.SIR.threshold.bootstrap = function(x,...){
    
    if(!inherits(x,"SIR.threshold.bootstrap")) stop("Only use with \"SIR.threshold.bootstrap\" obects")
    
    
    # Histogramme du nombre de variables sélectionnées par le modèle
    par(mfrow = (c(1, 2)))
    barplot((table(x$Nb.var.selec)/x$Nb.replications)*100, ylab = "percent", xlab = "Number of variables selected")
    title("Sizes of the reduced models")
    
    # Histogramme du nombre de fois ou chaque variable a été sélectionnée 
    barplot((x$effectif.var/x$Nb.replications)*100, names.arg = colnames(x$b.opt), ylab = "percent", xlab = "variable name")
    title("Variables selected in the reduced models")
    
}