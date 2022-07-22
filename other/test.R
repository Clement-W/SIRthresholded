library(SIRthresholded)
#https://cran.r-project.org/web/packages/MXM/vignettes/MMPC_tutorial.html
wine <-  read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data", header = FALSE) 
colnames(wine) <- c('Type', 'Alcohol', 'Malic', 'Ash', 'Alcalinity', 'Magnesium', 'Phenols', 'Flavanoids', 
                    'Nonflavanoids', 'Proanthocyanins', 'Color', 'Hue', 'Dilution', 'Proline')
wine <- wine[,-1] # Remove the class information (type of cultivars)

Y <- wine$Nonflavanoids
X <- wine[, -8] # remove target variable Nonflavanoids

res = SIR_threshold_opt(scale(Y),scale(X),H=5,n_lambda = 200)
res2 = SIR_threshold_bootstrap(scale(Y),scale(X),H=5,n_replications = 500,n_lambda=200) # à stocker quelque part
save(res2,file = "data/res2.RData")
