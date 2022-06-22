source("generationDonnees.R")
source("SIR.R")
source("cosCarre.R")

#====================================
# Méthode de simulation
#====================================

#N.simu = nombre de simulations 
#nn = nombre d'echantillon dans x
#pp = nombre de variables dans x
#pp1 = nombre de variables utiles dans x
#RRATIO = ratio var(eps)/var(Y))
simulation <- function(N.simu = 30, nn = 1000, pp = 50, pp1 = 10, RRATIO = 0.15, correlationX = FALSE, fct = "cubique", graphic = TRUE) {

    # contient les coscarré pour sir classique, soft puis hard dans les colonnes 
    # pour chaque simulation 
    mat.cos2.clas.hard.soft <- matrix(0, ncol = 4, nrow = N.simu)

    # Précision et recall pour les deux types de seuillage
    mat.precision.hard.soft <- matrix(0, ncol = 2, nrow = N.simu)
    mat.recall.hard.soft <- matrix(0, ncol = 2, nrow = N.simu)

    # variables sélectionnées pour le seuillage dur 
    liste.var.select.hard <- list()
    # variables sélectionnées pour le seuillage doux
    liste.var.select.soft <- list()
    cat("Starting ", N.simu, " simulations : ")

    lambdasHard = rep(0, N.simu)
    lambdasSoft = rep(0, N.simu)

    for (s in 1:N.simu) {
        #cat("s = ", s, " / ", N.simu)
        print(paste(s, "/", N.simu))
        # génération de données

        set.seed(s)
        don <- generationDonneesAvecRatio(n = nn, p = pp, p1 = pp1, graphic = graphic, RATIO = RRATIO, correlationX = correlationX, fct = fct)

        # Vecteur contenant les noms des dimensions utiles dans le modèle
        target.dims.to.keep = names(don$beta[which(don$beta != 0),])

        # Application de SIR classique
        res.SIR <- SIR.classique(don$Y, don$X, H = 10)
        # Applicatoin de SIR avec seuillage dur
        res.hard <- Sparse.SIR(don$Y, don$X, N.lambda = 200, graphic = FALSE, thresholding = "hard", output = FALSE)
        # Application de SIR avec seuillage doux
        res.soft <- Sparse.SIR(don$Y, don$X, N.lambda = 200, graphic = FALSE, thresholding = "soft", output = FALSE)

        # Application de SIR classique sur les variables sélectionnées par th soft
        res.final <- SIR.classique(don$Y, don$X[, res.soft$list.relevant.variables, drop = FALSE], H = 10)
        # on étend le beta pour qu'il ait une taille pp 
        beta.final = matrix(rep(0, pp), nrow = 1)
        colnames(beta.final) <- colnames(don$X)
        # On remplit les colonnes du beta.final avec les valeurs des colonnes du
        # beta estimé par SIR avec le nombre de variable réduit
        beta.final[which(colnames(beta.final) %in% colnames(res.final$beta))] = res.final$beta

        # Remplissage de la matrice de cos^2 pour évaluer l'estimation de beta
        mat.cos2.clas.hard.soft[s, 1] <- cosine_squared(res.SIR$beta, don$beta)
        mat.cos2.clas.hard.soft[s, 2] <- cosine_squared(res.hard$b.opt, don$beta)
        mat.cos2.clas.hard.soft[s, 3] <- cosine_squared(res.soft$b.opt, don$beta)
        mat.cos2.clas.hard.soft[s, 4] <- cosine_squared(beta.final, don$beta)

        # Remplissage de la liste des variables sélectionnées pour les 2 types de seuillage
        liste.var.select.hard[[s]] <- res.hard$list.relevant.variables
        liste.var.select.soft[[s]] <- res.soft$list.relevant.variables

        lambdasSoft[s] = res.soft$lambda.opt
        lambdasHard[s] = res.hard$lambda.opt


        # Remplissage de la matrice de Precision 
        # Precision = TP / TP+FP = TP / (nb de dimensions gardées par le modèle)
        TP.hard = sum(res.hard$list.relevant.variables %in% target.dims.to.keep)
        TP.soft = sum(res.soft$list.relevant.variables %in% target.dims.to.keep)

        mat.precision.hard.soft[s, 1] = TP.hard / length(res.hard$list.relevant.variables)
        mat.precision.hard.soft[s, 2] = TP.soft / length(res.soft$list.relevant.variables)

        # Remplissage de la matrice de Recall  
        # Recall = TP / TP+FN = TP / (nb de dimensions utiles à garder)
        mat.recall.hard.soft[s, 1] = TP.hard / length(target.dims.to.keep)
        mat.recall.hard.soft[s, 2] = TP.soft / length(target.dims.to.keep)

    }
    cat(fill = TRUE)

    if (graphic) {
        # Affichage des boxplots montrant la qualité de l'estimation des beta selon les méthodes
        par(mfrow = (c(2, 2)))
        res1 = boxplot(list(mat.cos2.clas.hard.soft[, 1], mat.cos2.clas.hard.soft[, 2], mat.cos2.clas.hard.soft[, 3], mat.cos2.clas.hard.soft[, 4]), names = c("SIR", "SIR with hard th", "SIR with soft th", "SIR on variables selected by soft th"))
        title(paste("cos^2 entre beta théorique et beta estimé selon les méthodes\n", "n=", nn, ", p=", pp, ", p1=", pp1, ", ratio=", RRATIO, "link fct=", fct, "CorrelationX=", correlationX))
        print("estimation coscarre :")
        print(res1$stats)
    }

    # vecteurs contenant le nombre de variables sélectionnées pour chaque simulation
    # pour les 2 méthodes de seuillage
    nb.var.select.hard <- rep(0, N.simu)
    nb.var.select.soft <- rep(0, N.simu)
    for (s in 1:N.simu) {
        nb.var.select.hard[s] <- length(liste.var.select.hard[[s]])
        nb.var.select.soft[s] <- length(liste.var.select.soft[[s]])
    }

    # Affichage de la sélection des variables
    if (graphic) {
        res2 = boxplot(nb.var.select.hard, nb.var.select.soft, names = c("hard", "soft"))
        title(paste("Nombre de variable sélectionnée sur ", N.simu, " simulations (avec", pp1, "variables utiles)"))
        print("selection de variables :")
        print(res2$stats)
    }


    # Boxplot de la précision et du recall
    if (graphic) {
        res3 = boxplot(mat.precision.hard.soft[, 1], mat.recall.hard.soft[, 1], mat.precision.hard.soft[, 2], mat.recall.hard.soft[, 2], names = c("precision hard", "recall hard", "precision soft", "recall soft"))
        title("Precision = Taux de dimensions gardées qui étaient réellement à sauvegarder \n Recall = Taux de dimensions à sauvegarder qui ont été réellement sauvegardées")
        print("boxplot precision recall :")
        print(res3$stats)
    }


    if (graphic) {
        plot(mat.recall.hard.soft[, 1], mat.precision.hard.soft[, 1], col = "blue", pch = 15, cex = 2, ylab = "precision", xlab = "recall", ylim = c(0, 1.4), xlim = c(0, 1))
        points(mat.recall.hard.soft[, 2], mat.precision.hard.soft[, 2], col = "orange", pch = 18, cex = 2)
        legend("topright", legend = c("Hard thresholding", "Soft thresholding"), col = c("blue", "orange"), pch = c(15, 18))
        title("Précision et recall")
    }

    return(list(precision = mat.precision.hard.soft, recall = mat.recall.hard.soft, betacos2 = mat.cos2.clas.hard.soft, listevarHard = liste.var.select.hard, listevarSoft = liste.var.select.soft, colnames = colnames(don$X), lambdasHard = lambdasHard, lambdasSoft = lambdasSoft))
}


#====================================
# Méthode de simulation
#====================================

nbSim <- 100 # nombre de simulations 
n <- 200 # nombre d'echantillon dans x
p <- 50 # nombre de variables dans x
p1 <- 10 # nombre de variables utiles dans x
RATIO <- 0.1 # ratio var(eps)/var(Y))
correlationX = FALSE
fct = "cubique"


a = simulation(N.simu = nbSim, nn = n, pp = p, pp1 = p1, RRATIO = RATIO, correlationX = correlationX, fct = fct, graphic = TRUE)
save(a, file = paste("simulation-nbsim", nbSim, "-n", n, "-p", p, "-p1:", p1, "-ratio", RATIO, ".RData", sep = ""))

# package pour output un tableau : knitr  


nbSim <- 100 # nombre de simulations 
correlationX = FALSE
fct = "cubique"
p.p1 = rbind(c(25, 50), c(5, 10))
for (i in 1:3) {
    p = p.p1[1, i]
    p1 = p.p1[2, i]
    for (varN in c(200, 300, 500)) {
        n = varN
        for (varRatio in c(0.1, 0.01)) {
            RATIO = varRatio

            print(paste("Simulation with p=", p, "; p1=", p1, "; n=", n, "; RATIO=", RATIO))
            resSimul = simulation(N.simu = nbSim, nn = n, pp = p, pp1 = p1, RRATIO = RATIO, correlationX = correlationX, fct = fct, graphic = TRUE)
            save(resSimul, file = paste("simulation-nbsim", nbSim, "-n", n, "-p", p, "-p1:", p1, "-ratio", RATIO, ".RData", sep = ""))
        }
    }
}

