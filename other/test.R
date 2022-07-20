
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
res = SIR_threshold_bootstrap(scale(target_NonFlav),scale(wine_dataset),H=5,n_replications = 500) # à stocker quelque part