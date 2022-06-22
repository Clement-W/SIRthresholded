source("generationDonnees.R")
source("SIR.R")
source("cosCarre.R")

#=================================
# Exemple
#------------
p1 = 10
p = 50
set.seed(71)
don <- generationDonneesAvecRatio(n = 200, p = p, p1 = p1, graphic = TRUE, RATIO = 0.1, correlationX = FALSE, fct = "cubique")

res.SIR <- SIR.classique(don$Y, don$X, H = 10)

b.SIR.bootstrap = SIR.bootstrap(don$Y, don$X, H = 10, B = 10)
res <- Sparse.SIR(don$Y, don$X, N.lambda = 200, graphic = TRUE, thresholding = "hard", beta.th = don$beta)


res2 <- Sparse.SIR.auto(don$Y, don$X, H = 10, Nb.replications = 100, graphic = TRUE,
                        N.lambda = 200, k = 2, thresholding = "hard")


res.SIR$beta
res$b.opt
res2$b.opt

res2$list.relevant.variables

# Modele simplifie avec p^* variables selectionnees
res.final <- SIR.classique(don$Y, don$X[, res2$list.relevant.variables], H = 10)
beta.final = matrix(rep(0, p), nrow = 1)
colnames(beta.final) <- colnames(don$X)
beta.final[which(colnames(beta.final) %in% colnames(res.final$b))] = res.final$beta

cosca.SIR.Classique = cosine_squared(res.SIR$b, don$beta)
cosca.Sparse.SIR = cosine_squared(res$b.opt, don$beta)
cosca.SIR.auto = cosine_squared(res2$b.opt, don$beta)
cosca.SIR.final = cosine_squared(beta.final, don$beta)

# Comparaison des directions estimees
cat("Cos^2 of estimated directions with target beta :", fill = TRUE)
cat(paste("classic SIR : ", cosca.SIR.Classique), fill = TRUE)
cat(paste("Sparse SIR : ", cosca.Sparse.SIR), fill = TRUE)
cat(paste("Sparse SIR auto (with multiple replications) : ", cosca.SIR.auto), fill = TRUE)
cat(paste("Final classic SIR on relevant variables : ", cosca.SIR.final), fill = TRUE)
