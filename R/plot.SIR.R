#' @export
#' @keywords internal
plot.SIR <- function(x, ...) {

    if (!inherits(x, "SIR"))
        stop("Only use with \"SIR\" obects")

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


    dev.new()
    plot(x$index_pred, x$Y, xlab = "Estimated first index", ylab = "y", pch = 4)
    title("Reconstructed index")

}