library(SIRthresholded)

a = read.csv("other/ozone.txt",sep = "")
Y = data.matrix(a[1],rownames.force = NA)
X = data.matrix(a[2:11],rownames.force = NA)

a=SIR_threshold_opt(Y,X,thresholding = "soft",n_lambda = 30)

b = read.csv("other/uscrime.txt",sep="\t")
Y = data.matrix(b[16],rownames.force = NA)
X = data.matrix(b[1:15],rownames.force = NA)

a=SIR_threshold_opt(scale(Y),scale(X),thresholding = "soft",n_lambda = 500)
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

res = SIR_threshold_opt(scale(target_NonFlav),scale(wine_dataset),H=10)


### https://journal.r-project.org/archive/2016/RJ-2016-009/RJ-2016-009.pdf
# FONCTIONNE BIEN:
load("other/pollution.RData")
X <- (pollution[, -19])
Y <- (pollution[, 19])

res = SIR_threshold_opt(Y,X,n_lambda = 300,thresholding = "soft")
