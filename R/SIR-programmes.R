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


#=============== Méthode SIR Classique =================================
# Y = Données y de R^1
# X = Données x de R^p
# H = nombre de tranche
SIR.classique <- function(Y, X, H = 10) {

    # Pour gérer le cas ou X ne contient qu'une variable
    if (is.null(dim(X)) || length(dim(X)) == 1) {
        X = matrix(X, ncol = 1)
    }

    n <- nrow(X)
    p <- ncol(X)


    # Moyenne des échantillons pour chaque variable de X
    moy.X <- matrix(apply(X, 2, mean), ncol = 1)

    # Calcul de la matrice de covariance de X
    Sigma <- var(X) * (n - 1) / n

    # Ordonne les données de X à partir des indices triés dans l'ordre croissant de Y
    X.ord <- X[order(Y),, drop = FALSE]

    # le vecteur nH contient le nombre d'observation pour chaque tranche H
    vect.nh <- rep(n %/% H, H) # %/% = integer division

    # si le nombre d'observations total est plus grand que le nb d'observation dans vect.nh
    # pour prendre le reste de la division de n %/% H 
    if ((n - sum(vect.nh)) > 0) {
        # pour chaque observation non présente dans vect.nh
        for (i in 1:(n - sum(vect.nh))) {
            # on prend une tranche au hasard
            h <- sample(1:H, 1)
            # on ajoute 1 observation dans cette tranche
            vect.nh[h] <- vect.nh[h] + 1
        }
    }

    # vect.ph = proportion des yi tombant dans chaque tranche
    vect.ph <- vect.nh / n

    # séparation des n échantillons en H tranche
    # (on assigne un indice h à chaque échantillon)
    vect.h <- rep(1:H, times = vect.nh)

    # Matrice des moyennes par tranche des xi (matrice p*H)
    mat.mh <- matrix(0, ncol = H, nrow = p)

    # Pout chaque tranche
    for (h in 1:H) {
        # récupération des indices pour la tranche h
        indicesTrancheH = which(vect.h == h)
        # récupération des x appartenant à la tranche h
        xDansTrancheH = X.ord[indicesTrancheH,, drop = FALSE]
        # applique la moyenne aux échantillons de chaque feature de la tranche h
        mat.mh[, h] <- apply(xDansTrancheH, 2, mean)
    }

    # Estimation de la matrice de covariances des moyennes par tranche
    moy.X.etendue = moy.X %*% matrix(rep(1, H), ncol = H) # Pour broadcaster l'opération - entre mat.mh et moy.X
    M <- (mat.mh - moy.X.etendue) %*% diag(vect.ph) %*% t(mat.mh - moy.X.etendue)

    # Matrice d'intérêt : inverse de la matrice de covariance de X multipliée par la
    # matrice de covariance des moyennes par tranche :(pxp) %*% (pxp)
    #  matrice d'intérêt =  Sigma^-1 * Cov(Moyenne par tranche)
    SigmaInv <- solve(Sigma)
    M1 <- SigmaInv %*% M

    # eigen(M1)$vectors donne une matrice pxp qui contient les p vecteurs propres de x dans chaque colonne
    # Les vecteurs propres sont retournés dans l'ordre décroissant de leur valeur propre associée.
    # Donc ici on récupère le vecteur propre associé à la plus grande valeur propre pour avoir b.est
    # donc la première colonne
    rSIR = eigen(M1)
    eig.values = Re(rSIR$values)
    eig.vectors = Re(rSIR$vectors)
    b.est <- eig.vectors[, 1]

    # conversion en matrice à une ligne et p colonnes
    b.est <- matrix(b.est, nrow = 1)

    # Ajout des noms de colonnes
    colnames(b.est) <- colnames(X)
    # On a donc b.est l'estimation de la direction de Beta

    res = list(beta = b.est, M1 = M1, eig.val = eig.values, eig.vect = eig.vectors, n = n, p = p)
    class(res) = "SIR"

    return(res)
}



#=========== Méthode SIR avec Bootstrap ============
# Y = Données y de R^1
# X = Données x de R^p
# H = nombre de tranche
# B = Nombre d'échantillon bootstrap à tirer
SIR.bootstrap <- function(Y, X, H = 10, B = 10) {
    p <- ncol(X)
    n <- nrow(X)

    # Initialisation de la matrice des estimations de beta_hat de taille pxB
    # avec l'estimation de beta_hat pour le nouvel échantillon après bootstrap dans
    # les colonnes de la matrice
    mat.b.est <- matrix(0, ncol = B, nrow = p)
    for (r in 1:B) {
        # indice des éléments à prendre dans X : on tire avec remise les indices des éléments de X
        # et on créé une liste qui contient ces indices (Bootstrap)
        indice <- sample(1:n, replace = TRUE)
        # on stock chaque estimation b dans les colonnes de la matrice d'estimation
        mat.b.est[, r] <- SIR.classique(Y[indice], X[indice,], H = H)$beta
    }

    # on récupère le vecteur propre associé à la plus grande valeur propre
    # de la matrice d'estimations de beta_hat
    b.boot <- eigen(mat.b.est %*% t(mat.b.est), symmetric = TRUE)$vectors[, 1]

    # conversion en matrice à une ligne
    b.boot <- matrix(b.boot, nrow = 1)

    # Renommage des colonnes
    colnames(b.boot) <- colnames(X)
    return(b.boot)
}



#==================== Méthode SIR avec seuillage avec un paramètre lambda ==================

# récupération des méthodes de seuillage :


# Y = Données y de R^1
# X = Données x de R^p
# H = nombre de tranche
# lambda = seuil pour le seuillage doux ou dur
# thresholding = méthode de seuillage utilisée
sparseSIR.lambda <- function(Y, X, H = 10, lambda = 0, thresholding = "hard") {
    n <- nrow(X)
    p <- ncol(X)

    # Estimation de la direction des beta et de la matrice d'intérêt Sigma^-1 * Cov(Moyenne par tranche)
    # avec la méthode SIR classique :
    res.clas <- SIR.classique(Y, X, H = H)
    b.clas <- res.clas$beta
    M1 <- res.clas$M1

    # Seuillage de la matrice d'intérêt avec la méthode indiquée en paramètre
    if (thresholding == "soft") {
        M1_th <- do_soft_thresholding(M1, lambda = lambda)
    }
    if (thresholding == "hard") {
        M1_th <- do_hard_thresholding(M1, lambda = lambda)
    }

    res.eig = eigen(M1_th)
    eig.values = Re(res.eig$values)
    eig.vectors = Re(res.eig$vectors)

    # Récupération du vecteur propre associé à la plus grande valeur propre
    # de la matrice d'intérêt seuillée
    b.estim <- eig.vectors[, 1]

    # nombre de zéros présents dans l'estimation du b
    nb.zeros <- length(which(b.estim == 0))

    # Calcul de la qualité de la corrélation entre les b estimés par la méthode classqiue
    # et par le sir thresholding
    cosca <- cosine.squared(b.clas, b.estim)

    # Conversion en matrice en une ligne
    b.estim <- matrix(b.estim, nrow = 1)
    colnames(b.estim) <- colnames(X)

    # Liste des variables utiles qui sont les colonnes de l'estimation de b
    # qui sont différentes de 0
    if (nb.zeros == 0) {
        list.relevant.variables = colnames(X)
    } else {
        list.relevant.variables <- colnames(b.estim)[-which(b.estim == 0)]
    }


    res = list(beta = b.estim, M1_th = M1_th, eig.val = eig.values, eig.vect = eig.vectors, n = n, p = p, nb.zeros = nb.zeros, list.relevant.variables = list.relevant.variables, cosca = cosca)
    class(res) = "SIR-lambda"

    return(res)

}


#=========== Méthode SIR avec seuillage et N.lambda paramètres lambda =========


# Y = Données y de R^1
# X = Données x de R^p
# H = nombre de tranche
# N.lambda = nombre de lambda à tester
# thresholding = méthode de seuillage utilisée
Sparse.SIR <- function(Y, X, H = 10, N.lambda = 100, thresholding = "hard", graphic = TRUE, output = TRUE, beta.th = NULL) {

    n <- nrow(X)
    p <- ncol(X)

    # Estimation de la direction des beta et de la matrice d'intérêt  M1 = Sigma^-1 * Cov(Moyenne par tranche)
    # avec la méthode SIR classique :
    resSIR <- SIR.classique(Y, X, H = 10)
    b.clas <- resSIR$beta
    M1 <- resSIR$M1

    # Création d'une liste de lambdas allant de 0 à la valeur absolue maximum de M1, avec un total
    # de 100 valeurs.
    lambdas <- seq(0, max(abs(M1)), length.out = N.lambda + 1)[-(N.lambda + 1)]

    # Initialisation d'une matrice de 0 de taille N.lambda*p
    # pour reccueillir les estimations du vecteur b pour chacun des lambdas
    mat.b.th <- matrix(0, ncol = p, nrow = N.lambda)


    # Initialisation d'un vecteur de N.lambda NA pour receuillir le nombre de 0 après seuillage
    # pour chaque valeur de lambda dans l'estimation du vecteur beta
    vect.nb.zeros <- rep(NA, N.lambda)

    # Initialisation d'un vecteur de N.lambda NA pour recueillir les cos carrés entre b et b_seuillage
    # pour chaque valeur de lambda
    vect.cosca <- rep(NA, N.lambda)
    vect.cosca.final <- rep(NA, N.lambda)

    if (!is.null(beta.th)) {
        vect.cosca.th = rep(NA, N.lambda)
        vect.cosca.final.th = rep(NA, N.lambda)
    }

    # Initialisation de la liste des variables utiles
    list.relevant.variables <- list()

    b.p1 = list()

    # Application de la méthode SIR avec seuillage, avec les N.lambda valeurs
    # de lambdas
    for (i in 1:N.lambda) {

        # Récupération du résultat après application de la méthode SIR avec le lambda_i
        resSparseSIR <- sparseSIR.lambda(Y, X, H = H, lambda = lambdas[i], thresholding = thresholding)

        # Stockage de l'estimation du vecteur beta dans la ligne i de la matrice
        mat.b.th[i,] <- resSparseSIR$beta

        # Stockage du nombre de 0 dans l'estimation du vecteur beta après seuillage pour
        # cette valeur de lambda
        vect.nb.zeros[i] <- resSparseSIR$nb.zeros

        # Stockage du cos carré entre b et b_seuillage pour cette valeur de lambda
        vect.cosca[i] <- resSparseSIR$cosca


        if (!is.null(beta.th)) {
            vect.cosca.th[i] = cosine.squared(beta.th, resSparseSIR$beta)
        }

        if (length(resSparseSIR$list.relevant.variables) != 0 && !is.null(beta.th)) {
            res.final <- SIR.classique(Y, X[, resSparseSIR$list.relevant.variables, drop = FALSE], H = 10)


            # on étend le beta pour qu'il ait une taille pp 
            beta.final = matrix(rep(0, p), nrow = 1)
            colnames(beta.final) <- colnames(X)
            # On remplit les colonnes du beta.final avec les valeurs des colonnes du
            # beta estimé par SIR avec le nombre de variable réduit
            beta.final[which(colnames(beta.final) %in% colnames(res.final$b))] = res.final$b

            vect.cosca.final.th[i] = cosine.squared(beta.th, beta.final)
            vect.cosca.final[i] = cosine.squared(b.clas, beta.final)

            b.p1[[i]] = res.final$b
        }


        # Stocage des variables utiles pour cette valeur de lambda à l'indice i de la liste
        # des variables utiles
        list.relevant.variables[[i]] <- resSparseSIR$list.relevant.variables
    }



    # Création d'un vecteur qui contient p valeurs : pour chaque variable x_i
    # est associé l'indice du lambda à partir duquel la variable devient inutile. 
    # Cet indice correspond donc également au nombre de lambda pour lesquels la variable
    # est utile.
    indice.0 = colSums(mat.b.th / mat.b.th, na.rm = TRUE)

    # On recherche ensuite un point de rupture dans la liste indice.0 ordonnée. Ce point de rupture 
    # fixe le nombre de variables à ne pas sélectionner, soit le nombre de variable 
    # où l'estimation de beta donne zero pour un lambda donné.
    require(strucchange)
    fit_bp <- breakpoints(sort(indice.0, decreasing = FALSE) ~ 1, breaks = 1, h = 2 / p)

    # Si le nombre de coeficients de l'estimation de beta égaux à 0 et qui sont égaux au point
    # de rupture, donc au nombre de variable à ne pas sélectionner : 
    if (length(which(vect.nb.zeros == fit_bp$breakpoints)) > 0) {
        # l'indice optimal de lambda est l'indice minimum où le nombre de 0 apparaissant pour une valeur
        # de lambda est égal au nombre de valeurs inutiles au point de rupture
        indice.opt <- min(which(vect.nb.zeros == fit_bp$breakpoints))
    } else {
        # Sinon, on somme le nombre de fois où le nombre de 0 qui apparait dans l'estimation
        # de b pour chaque lambda est inférieur au nombre de valeurs inutiles au point de rupture
        # + 1
        indice.opt <- sum(vect.nb.zeros < fit_bp$breakpoints) + 1
    }

    # On récupère le lambda qui correspond à cet indice
    lambda.opt <- lambdas[indice.opt]

    # Affichage de la sélection du lambda
    if (graphic == TRUE) {

        if (length(which(vect.nb.zeros == fit_bp$breakpoints)) > 0) {
            #par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
            # affichage des indices triés
            plot(sort(indice.0, decreasing = FALSE), xlab = "variable i", ylab = expression("Indices of" ~ lambda))
            # ligne verticale pour montrer la rupture
            abline(v = fit_bp$breakpoints + 0.5, col = 6, lwd = 3)
            title(paste("Choix de l'indice du lambda optimal parmi les", N.lambda, "testés"))
        }
        else {
            par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
            # affichage des indices triés
            plot(sort(indice.0, decreasing = FALSE), xlab = "variable p_i", ylab = "Indice du lambda à partir duquel la variable est inutile")
            # ligne verticale pour montrer la rupture
            abline(v = fit_bp$breakpoints + 1.5, col = 7, lwd = 3)
            title(paste("Choix de l'indice du lambda optimal parmi les ", N.lambda, " testés"))
        }

        # Affichage du pourcentage de variable utiles en fonction des lambdas
        #     plot(lambdas, 1 - vect.nb.zeros / p, ylim = c(0, 1.3), xlim = c(0, max(lambdas)),xlab = expression(lambda), col = 3, type = "l", ylab = "",lwd=3)
        # 
        #     # Affichage de l'évolution du cos carré entre SIR et SParseSIR en fonction des valeurs de lambda
        #     points(lambdas, vect.cosca, type = "l", ylab = "", xlab = "lambda",lwd=3,col="red")
        #     title(paste("Evolution of cos^2 as a function of lambda (with",thresholding,"thresholding)"))
        #     #print(vect.cosca[indice.opt])
        # 
        #     if (!is.null(beta.th)) {
        #         # Affichage de l'évolution du cos2 entre SparseSIR et betaTh
        #         #points(lambdas, vect.cosca.th, type = "l", ylab = " ", ylim = c(0, 1.2), xlab = "lambda", col = "red")
        # 
        #         # pareil avec final
        #         #points(lambdas, vect.cosca.final.th, type="l", ylab = " ", ylim = c(0, 1.2), xlab = "lambda", col = "blue")
        #         
        #         # Affichage de l'évolution du cos2 entre SIR et beta.final
        #         # COURBE JAUNE À INCLURE À LA FIN
        #         #points(lambdas, vect.cosca.final, type = "l", ylab = " ", ylim = c(0, 1.2), xlab = "lambda", col = "orange",lwd=3)
        #         
        #         
        #         # Ligne verticale au maximum du cos carré entre le beta théorique et le beta estimé par SparseSIR
        #         abline(v = lambdas[which(vect.cosca.th == max(vect.cosca.th))], col = "orange", lty = 3,lwd=3)
        # 
        #         # Ligne horizontale au cos carré entre le beta théorique et le beta estimé par SIR classique
        #         abline(h = cosine.squared(beta.th, b.clas), col = "gray60", lty = 3,lwd=3)
        #         #print(cosine.squared(beta.th, b.clas))
        #         
        #         #print(cosine.squared(b.p1[[indice.opt]],rep(1,p1)))
        #         #print(vect.cosca.final.th[indice.opt])
        #        
        #     }
        #     # Ligne verticale au lambda optimal
        #     abline(v = lambda.opt, col = 6, lwd = 3)
        # 
        #     # Ajout de la légende
        #     if (is.null(beta.th)) {
        #         legend("topright", legend = c("cos2(SIR,SparseSIR)", "cos2(BetaTh,SparseSIR)", "% of relevant var.", "optimal lambda"), col = c(1, "black", 3, 6), lty = 1)
        #     }
        #     else {
        #         #legend("topright", legend = c("cos2(SIR,SparseSIR)", "cos2(BetaTh,SparseSIR)","cos2(BetaTh,th-p*->p)","cos2(SIR,th-p*->p)", "% of relevant var.", "optimal lambda", "cos2(BetaTh,SIR)", "max(cos2(BetaTh,SparseSIR))"), col = c(1, "red", "blue","orange",3, 6, "gray60", "red"), lty = c(1, 1,1, 1,1, 1, 3, 3))
        #         
        #         #Part 1
        #         #legend("topright", legend = c("% of relevant variables", expression("optimal"~lambda), expression(cos^2~(beta~","~hat(b)[SIR]))), col = c(3, 6, "gray60"), lty = c(1, 1, 3),lwd=c(3,3,3),cex=1.6)
        #         
        #         #Part 2
        #         #legend("topright", legend = c("% of relevant variables", expression("optimal"~lambda), expression(cos^2~(beta~","~hat(b)[SIR])), expression(max~of~cos^2~(beta~","~hat(b)[HT-SIR])),expression(cos^2~(hat(b)[SIR]~","~hat(b)[HT-SIR]))), col = c(3, 6, "gray60", "blue","black"), lty = c( 1, 1, 3,3,1),lwd=c(3,3,3,3,3),cex=0.75)
        # 
        #         #Part 3
        #         #legend("topright", legend = c("% of relevant variables", "Optimal lambda", "cos²(trueBeta,beta_SIR)", "maximum of cos²(trueBeta,beta_HT-SIR)","cos²(beta_SIR,beta_HT-SIR)","cos²(beta_SIR,beta_HT-SIR-p*)"), col = c(3, 6, "gray60", "blue",1,"orange"), lty = c( 1, 1, 3,3,1,1),lwd=c(3,3,3,3,3,3),cex=0.75)
        #         
        #         
        #         #Comparison
        #         legend("topright", legend = c(expression(max~of~cos^2~(beta~","~hat(b)[ST-SIR])),expression(cos^2~(hat(b)[SIR]~","~hat(b)[ST-SIR]))), col = c("orange","red"), lty = c(3,1),lwd=c(3,3),cex=1.4)
        #         
        #     }
    }

    if (is.na(lambda.opt) == TRUE) {
        b.opt <- SIR.classique(Y, X, H = 10)$beta
        list.relevant.var <- colnames(X)
        if (output == TRUE) {
            affichageResultat(sparseSIR = FALSE, b.opt = b.opt, p = p)
        }
    }
    else {
        # Cas où on trouve bien un lambda optimal parmi les N.lambda
        # Si le nombre de variables inutiles est inférieur au nombre de variable total-1
        if (fit_bp$breakpoints < (p - 1)) {

            # Le beta optimal est la ligne de la matrice d'estimations de beta qui 
            # correspond à lambda.opt
            b.opt <- mat.b.th[which(lambdas == lambda.opt),]
            # Conversion du beta en matrice à une ligne p colonnes
            b.opt <- matrix(b.opt, nrow = 1)
            # Renommage des colonnes
            colnames(b.opt) <- colnames(X)
            # Récupération des variables utiles qui sont les colonnes du b.opt
            # où la valeur est différente de zero
            list.relevant.var <- colnames(b.opt)[-which(b.opt == 0)]

            if (output == TRUE) {
                affichageResultat(b.opt = b.opt, list.relevant.var = list.relevant.var, p = p)
            }
        }
        # Cas où la méthode Sparse SIR avec seuillage n'a pas permis de réaliser de la sélection de variable
        else {
            b.opt <- SIR.classique(Y, X, H = 10)
            if (output == TRUE) {
                affichageResultat(sparseSIR = FALSE, b.opt = b.opt, p = p)
            }
        }
    }


    res = list(b.opt = b.opt, lambdas = lambdas, lambda.opt = lambda.opt, mat.b.th = mat.b.th,
       N.lambda = N.lambda, vect.nb.zeros = vect.nb.zeros,
       list.relevant.variables = list.relevant.var)

    class(res) = "SIR-lambdas"
    return(res)
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
