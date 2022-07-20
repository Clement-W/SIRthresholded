library(SIRthresholded)

a = read.csv("other/ozone.txt",sep = "")
Y = data.matrix(a[1],rownames.force = NA)
X = data.matrix(a[2:11],rownames.force = NA)

a=SIR_threshold_opt(Y,X,thresholding = "soft",n_lambda = 30)

b = read.csv("other/uscrime.txt",sep="\t")
Y = data.matrix(b[16],rownames.force = NA)
X = data.matrix(b[1:15],rownames.force = NA)

a=SIR_threshold_opt(Y,X,thresholding = "soft",n_lambda = 500)
b = SIR(Y,X)

#https://cran.r-project.org/web/packages/MXM/vignettes/MMPC_tutorial.html
wine <-  read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data", header = FALSE) 
colnames(wine) <- c('Type', 'Alcohol', 'Malic', 'Ash', 
                    'Alcalinity', 'Magnesium', 'Phenols', 
                    'Flavanoids', 'Nonflavanoids', 'Proanthocyanins',
                    'Color', 'Hue', 'Dilution', 'Proline')

# 3. REMOVE the 1st attribute, which is the class information:
wine <- wine[,-1] 

targetVariable <- wine$Nonflavanoids
target_NonFlav <- as.matrix(wine$Nonflavanoids,ncol=1)
wine_dataset <- as.matrix(wine[, -8]) # remove target variable
wine_dataset[, 12] <- as.numeric(wine_dataset[, 12]) # change type int to num

res = SIR_threshold_opt(scale(target_NonFlav),scale(wine_dataset),H=5,n_lambda = 200)


### https://journal.r-project.org/archive/2016/RJ-2016-009/RJ-2016-009.pdf
# FONCTIONNE BIEN:
load("other/pollution.RData")
X <- (pollution[, -19])
Y <- (pollution[, 19])

res = SIR_threshold_opt(Y,X,n_lambda = 300,thresholding = "soft")


simulation <- function(N.simu = 30, nn = 1000, pp = 50, pp1 = 10, RRATIO = 0.15, correlationX = FALSE, fct = "cubique", graphic = TRUE) {
    
    # contient les coscarré pour sir classique, soft puis hard dans les colonnes 
    # pour chaque simulation 
    mat.cos2.clas.hard.soft <- matrix(0, ncol = 5, nrow = N.simu)
    
    # Précision et recall pour les deux types de seuillage
    mat.precision.hard.soft <- matrix(0, ncol = 3, nrow = N.simu)
    mat.recall.hard.soft <- matrix(0, ncol = 3, nrow = N.simu)
    
    # variables sélectionnées pour le seuillage dur 
    liste.var.select.hard <- list()
    # variables sélectionnées pour le seuillage doux
    liste.var.select.soft <- list()
    liste.var.select.scad = list()
    cat("Starting ", N.simu, " simulations :\n")
    
    lambdasHard = rep(0,N.simu)
    lambdasSoft = rep(0,N.simu)
    lambdasScad = rep(0,N.simu)
    
    # Réalisation des simulations
    for (s in 1:N.simu) {
        
        if (s%%5 == 0) {
            print(paste("simul", s, "/", N.simu))
        }
        
        # génération de données
        set.seed(s)
        don <- generationDonneesAvecRatio(n = nn, p = pp, p1 = pp1, graphic = graphic, RATIO = RRATIO, correlationX = correlationX, fct = fct)
        
        # Vecteur contenant les noms des dimensions utiles dans le modèle
        target.dims.to.keep = names(don$beta[which(don$beta != 0),])
        
        # Application de SIR classique
        res.SIR <- SIR(don$Y, don$X, H = 10,graph = FALSE)
        # Applicatoin de SIR avec seuillage dur
        res.hard <- SIR_threshold_opt(don$Y, don$X,n_lambda = 200, graph = FALSE, thresholding = "hard", output = FALSE)
        # Application de SIR avec seuillage doux
        res.soft <- SIR_threshold_opt(don$Y, don$X,n_lambda = 200, graph = FALSE, thresholding = "soft", output = FALSE)
        
        res.scad <- SIR_threshold_opt(don$Y, don$X,n_lambda = 200, graph = FALSE, thresholding = "scad", output = FALSE)
        
        # Application de SIR classique sur les variables sélectionnées par th soft
        res.final <- SIR(don$Y, don$X[, res.soft$list_relevant_variables, drop = FALSE], H = 10,graph = FALSE)
        # on étend le beta pour qu'il ait une taille pp 
        beta.final = matrix(rep(0, pp), nrow = 1)
        colnames(beta.final) <- colnames(don$X)
        # On remplit les colonnes du beta.final avec les valeurs des colonnes du
        # beta estimé par SIR avec le nombre de variable réduit
        beta.final[which(colnames(beta.final) %in% colnames(res.final$b))] = res.final$b
        
        # Remplissage de la matrice de cos^2 pour évaluer l'estimation de beta
        mat.cos2.clas.hard.soft[s, 1] <- cosine_squared(res.scad$b, don$beta)
        mat.cos2.clas.hard.soft[s, 2] <- cosine_squared(res.hard$b, don$beta)
        mat.cos2.clas.hard.soft[s, 3] <- cosine_squared(res.soft$b, don$beta)
        mat.cos2.clas.hard.soft[s, 4] <- cosine_squared(beta.final, don$beta)
        
        # Remplissage de la liste des variables sélectionnées pour les 2 types de seuillage
        liste.var.select.hard[[s]] <- res.hard$list_relevant_variables
        liste.var.select.soft[[s]] <- res.soft$list_relevant_variables
        liste.var.select.scad[[s]] <- res.scad$list_relevant_variables
        
        lambdasSoft[s] = res.soft$lambda_opt
        lambdasHard[s] = res.hard$lambda_opt
        lambdasScad[s] = res.scad$lambda_opt
        
        
        # Remplissage de la matrice de Precision 
        # Precision = TP / TP+FP = TP / (nb de dimensions gardées par le modèle)
        TP.hard = sum(res.hard$list_relevant_variables %in% target.dims.to.keep)
        TP.soft = sum(res.soft$list_relevant_variables %in% target.dims.to.keep)
        TP.scad = sum(res.scad$list_relevant_variables %in% target.dims.to.keep)
        
        mat.precision.hard.soft[s, 1] = TP.hard / length(res.hard$list_relevant_variables)
        mat.precision.hard.soft[s, 2] = TP.soft / length(res.soft$list_relevant_variables)
        mat.precision.hard.soft[s, 3] = TP.scad / length(res.scad$list_relevant_variables)
        
        # Remplissage de la matrice de Recall  
        # Recall = TP / TP+FN = TP / (nb de dimensions utiles à garder)
        mat.recall.hard.soft[s, 1] = TP.hard / length(target.dims.to.keep)
        mat.recall.hard.soft[s, 2] = TP.soft / length(target.dims.to.keep)
        mat.recall.hard.soft[s, 3] = TP.scad / length(target.dims.to.keep)
        
    }
    cat(fill = TRUE)
    
    
    # Affichage des sorties graphiques
    if (graphic) {
        
        # Affichage des boxplots montrant la qualité de l'estimation des beta selon les méthodes
        par(mfrow = (c(2, 3)))
        boxplot_cos2(betacos2 = mat.cos2.clas.hard.soft,title=expression(cos ^ 2 ~ (hat(b) ~"," ~ beta)))
        
        # Affichage de la sélection des variables
        nbVarSelectionnees(listevarHard = liste.var.select.hard, listevarSoft = liste.var.select.scad,title="Number of selected variables")
        
        # Boxplot de la précision 
        boxplotPrecision(precision = mat.precision.hard.soft,title="Precision")
        
        #Boxplot du recall
        boxplotRecall(recall = mat.recall.hard.soft,title="Recall")
        
        #Boxplot occurrences
        barplot_occurrences(listevarHard = liste.var.select.hard,listevarSoft = liste.var.select.soft,title="% of variable selection")
        
        #Boxplot model size
        barplot_model_size(listevarHard = liste.var.select.hard,listevarSoft = liste.var.select.scad,"Model size")

    } 
    
    return(list(precision = mat.precision.hard.soft, recall = mat.recall.hard.soft, betacos2 = mat.cos2.clas.hard.soft, listevarHard = liste.var.select.hard, listevarSoft = liste.var.select.soft, colnames = colnames(don$X), lambdasHard=lambdasHard, lambdasSoft=lambdasSoft))
}

nbSim <- 20 # nombre de simulations 
n <- 200 # nombre d'echantillon dans x
p <- 50 # nombre de variables dans x
p1 <- 10 # nombre de variables utiles dans x
RATIO <- 0.1 # ratio var(eps)/var(Y))
correlationX = FALSE
fct = "cubique"

resSimul = simulation(N.simu = nbSim, nn = n, pp = p, pp1 = p1, RRATIO = RATIO, correlationX = correlationX, fct = fct, graphic = TRUE)

