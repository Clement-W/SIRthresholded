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









#=========== Méthode sparse SIR sur N réplications =========
# Effectue la méthode sparse SIR avec thrésholding sur N réplications

#k = facteur de multiplication de la taille de l'échantillon bootstrapé
Sparse.SIR.auto <- function(Y, X, H = 10, thresholding = "hard", Nb.replications = 200, graphic = TRUE, output = TRUE, N.lambda = 50, k = 1) {

    # Sparse SIR avec N.lambda sur tout l'échantillon
    res.SparseSIR <- Sparse.SIR(Y, X, H = H, thresholding = thresholding, graph = FALSE, output = FALSE, N.lambda = N.lambda)

    p <- ncol(X)
    n <- nrow(X)

    # Stockage du nombre de variable sélectionnée pour chaque réplication
    Nb.var.selec <- rep(NA, Nb.replications)
    # vecteur des noms de variables de X
    list.relevant.variables <- colnames(X)

    # Tableau permettant de stocker les variables sélectionnées pour chaque réplication
    liste <- list()

    # Vecteur qui contient le nombre de fois ou la variable à l'indice j a été retenue
    # dans une réplication
    effectif.var <- rep(0, p)

    # Pour chaque réplication
    for (replic in 1:Nb.replications) {
        # bootstrap 
        indice.boot <- sample(1:n, size = k * n, replace = TRUE)
        # On récupère les X et Y avec les indices pris par bootstrap
        X.boot <- X[indice.boot,]
        Y.boot <- Y[indice.boot]

        # Sparse.SIR sur l'echantillon créé par boostrap
        res.boot <- Sparse.SIR(Y.boot, X.boot, H = H, thresholding = thresholding, graph = FALSE, output = FALSE, N.lambda = N.lambda)
        # Stockage du nombre de variable sélectionnée (utiles) pour cette réplication
        Nb.var.selec[replic] <- length(res.boot$list.relevant.variables)
        print(length(res.boot$list.relevant.variables))

        liste[[replic]] <- res.boot$list.relevant.variables # stockage des variables sélectionnées
    }

    # Pour chaque réplication
    for (i in 1:Nb.replications) {
        # Pour chaque variable
        for (j in 1:p) {
            # Si la variable à l'indice j a été retenue pour la ième réplication
            if ((sum(liste[[i]] == colnames(X)[j])) == 1) {
                # On incrémente le compteur d'occurrence de la variable j
                effectif.var[j] <- effectif.var[j] + 1
            }
        }
    }

    if (graphic == TRUE) {
        # Histogramme du nombre de variables sélectionnées par le modèle
        par(mfrow = (c(1, 2)))
        barplot(table(Nb.var.selec), ylab = "Number of replication", xlab = "Number of variables selected ")

        # Histogramme du nombre de fois ou chaque variable a été sélectionnée 
        barplot(effectif.var, names.arg = colnames(X), ylab = "Number of replication where the variable has been selected", xlab = "variable name")
    }

    # Nombre de variable sélectionnée optimale (la variable qui a été la plus sélectionnée au total 
    # sur toutes les réplications effectuées)
    Nb.var.selec.opt <- as.numeric(names(which(table(Nb.var.selec) == max(table(Nb.var.selec)))))

    # Nombre de zero optimal
    Nb.zeros.opt <- p - Nb.var.selec.opt

    # estimation du beta final en prenant le beta estimé sur tout l'échantillon par la
    # méthode SIR, au lambda à partir duquel le nombre de zero optimal apparaît
    b.opt.final <- res.SparseSIR$mat.b.th[min(which(res.SparseSIR$vect.nb.zeros == Nb.zeros.opt)),]
    # Conversion du beta en matrice à une ligne p colonnes
    b.opt.final <- matrix(b.opt.final, nrow = 1)
    # Renommage des colonnes
    colnames(b.opt.final) <- colnames(X)

    # Si on a bien réduit le nombre de variables
    print(Nb.zeros.opt)
    if (Nb.zeros.opt > 0) {
        # mise à jour des variables utiles 
        list.relevant.variables <- colnames(X)[-which(b.opt.final == 0)]
        # Si le nombre de zéros dans le b.opt final est à 0
        if (length(which(b.opt.final == 0)) == 0) {
            # la liste de variables utile contient toutes les variables
            list.relevant.variables <- colnames(X)
        }
    }

    lambda.optim = res.SparseSIR$lambdas[min(which(res.SparseSIR$vect.nb.zeros == Nb.zeros.opt))]

    affichageResultat(b.opt = b.opt.final, list.relevant.var = list.relevant.variables, p = p)

    return(list(b.opt = b.opt.final, lamba.opt = lambda.optim,
       Nb.var.selec.opt = Nb.var.selec.opt, list.relevant.variables = list.relevant.variables))
}
