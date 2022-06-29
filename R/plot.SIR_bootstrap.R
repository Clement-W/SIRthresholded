#' @export
#' @keywords internal
plot.SIR_bootstrap <- function(x, choice="",...) {
    
    if (!inherits(x, "SIR_bootstrap"))
        stop("Only use with \"SIR_bootstrap\" obects")
    
    if (!(choice %in% c("eigvals", "estim_ind",""))) 
        stop("\"choice\" must be either \"eigvals\" or \"estim_ind\"",call. = FALSE)
    
    if(choice=="" || choice=="eigvals"){
        dev.new()
        # on prend au plus les 10 premières valeurs propres et on les affichent 
        # dans un screplot
        eig_vals <- x$eig_val
        title = "Eigenvalues"
        if (length(eig_vals) >= 10) {
            eig_vals <- eig_vals[1:10]
            title = "10 first eigenvalues"
        }
        
        a=barplot(eig_vals, ylab = "eigenvalues", xlab = "dimensions", main = title)
        axis(1, at = a, labels=seq(1,length(eig_vals))) 
    }
    
    
    if(choice=="" || choice=="estim_ind"){
        dev.new()
        plot(x$index_pred, x$Y, xlab = "Estimated first index", ylab = "y", pch = 4)
        title("Reconstructed index")
    }
    
}