#========== Affichage des résultats après méthod SIR)

affichageResultat <- function(sparseSIR = TRUE, b.opt, list.relevant.var = NULL, p) {
    cat("===============================", fill = TRUE)
    cat(paste("Estimated EDR direction by ", ifelse(sparseSIR, "SparseSIR :", "usual SIR (SparseSIR not possible!)")), fill = TRUE)
    cat(paste("EDR direction has been estimated by ", ifelse(sparseSIR, "SparseSIR :", "usual SIR :")), fill = TRUE)
    print(round(b.opt, digits = 5))
    cat("===============================", fill = TRUE)

    nb.relevant.var = ifelse(!is.null(list.relevant.var), length(list.relevant.var), p)

    cat(paste("Number of relevant variables = ", nb.relevant.var, " over the ", p, " available variables"), fill = TRUE)
    cat("===============================", fill = TRUE)

    if (!is.null(list.relevant.var)) {
        cat(paste("List of relevant variables:", paste(list.relevant.var, collapse = ",")), fill = TRUE)
    }
    cat("===============================", fill = TRUE)
}

