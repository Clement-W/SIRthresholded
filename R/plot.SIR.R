#' @export
#' @keywords internal
plot.SIR <- function(x, ...) {

    if (!inherits(x, "SIR"))
        stop("Only use with \"SIR\" obects")

    dev.new()
    # on prend au plus les 10 premières valeurs propres et on les affichent 
    # dans un screplot
    eig_vals <- x$eig_val
    title = "Eigen values"
    if (length(eig_vals) >= 10) {
        eig_vals <- eig_vals[1:10]
        title = "10 first eigen values"
    }
    plot(eig_vals, ylab = "eigen values", xlab = "", type = "l", main = title)
    points(eig_vals, pch = 16)

    dev.new()
    plot(x$index_pred, x$Y, xlab = "Estimated index", ylab = "y", pch = 4)
    title("Reconstructed index")



}