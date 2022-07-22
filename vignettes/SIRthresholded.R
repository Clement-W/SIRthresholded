## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo=TRUE,eval=TRUE,comment="#>")
custom_lambda = 0.071


## -----------------------------------------------------------------------------
set.seed(4)
n <- 200 # Sample size
p <- 30 # Number of variables in X
p_star <- 10 # Number of relevant variables in X
X <- mvtnorm::rmvnorm(n,sigma=diag(p)) # X ~ N(0,I_p)
dimnames(X) <- list(1:n, paste("X", 1:p, sep = "")) # Rename columns of X

eps <- rnorm(n, sd = 20) # Error

beta <- matrix(c(rep(1,p_star),rep(0,p-p_star)),ncol=1) # Beta = Heaviside function
rownames(beta) <- colnames(X) # Rename rows of X

Y <- (X %*% beta)**3 + eps # The model

## ----fig.align='center',fig.width=5,fig.asp=1,out.width="300px"---------------
par(mar=c(5,4,1,1)+0.1)
plot(X %*% beta, Y, xlab = "true index")

## -----------------------------------------------------------------------------
library(SIRthresholded)

## -----------------------------------------------------------------------------
res_SIR = SIR(Y = Y, X = X, H = 10,graph = FALSE)

## -----------------------------------------------------------------------------
summary(res_SIR)

## -----------------------------------------------------------------------------
cor(c(beta),c(res_SIR$b))

## ----fig.asp=1,fig.show='hold',fig.width=5,out.width="300px",fig.align='center',fig.cap=' '----
plot(res_SIR,choice="estim_ind")
plot(res_SIR,choice="eigvals")

## -----------------------------------------------------------------------------
res_STSIR = SIR_threshold_opt(Y=Y, X=X, H=10, n_lambda=100, thresholding="soft", graph=FALSE)

## -----------------------------------------------------------------------------
summary(res_STSIR)

## ----fig.width=7,out.width="450px",fig.asp=1,fig.align='center'---------------
plot(res_STSIR,choice="cos2_selec")

## ----fig.width=7,out.width="450px",fig.asp=1,fig.align='center'---------------
plot(res_STSIR,choice="opt_lambda")

## ----fig.width=7,out.width="500px",fig.asp=1,fig.align='center'---------------
plot(res_STSIR,choice="regul_path")

## ----fig.align='center',fig.width=5,fig.asp=1,out.width="300px"---------------
plot(res_STSIR,choice="estim_ind")

## ----fig.asp=1,fig.show='hold',fig.width=8,out.width="330px",fig.align='center',fig.cap=' '----
res_HTSIR = SIR_threshold_opt(Y=Y, X=X, H=10, n_lambda=100, thresholding="hard", graph=FALSE)

plot(res_HTSIR,choice="cos2_selec")
plot(res_HTSIR,choice="regul_path")

## -----------------------------------------------------------------------------
res_SIR_pstar = SIR(Y=res_STSIR$Y, X=res_STSIR$X_reduced, H=res_STSIR$H)
summary(res_SIR_pstar)

## ----fig.align='center',fig.width=5,fig.asp=1,out.width="300px"---------------
plot(res_SIR_pstar,choice="estim_ind")

## -----------------------------------------------------------------------------
b_extended <- matrix(rep(0,p),nrow=1) # Create the empty vector
colnames(b_extended) <- colnames(X) # Assign the colnames of X
# Assign the values of b_extended according to the colnames
b_extended[which(colnames(b_extended) %in% colnames(res_SIR_pstar$b))] = res_SIR_pstar$b

## -----------------------------------------------------------------------------
cor(c(beta),c(b_extended))

## -----------------------------------------------------------------------------
res_SIR_thresh = SIR_threshold(Y, X, H = 10, lambda = 0.04, thresholding = "hard")
summary(res_SIR_thresh)

## -----------------------------------------------------------------------------
res_SIR_boot = SIR_threshold_bootstrap(Y,X,H=10,n_lambda=100,thresholding="hard", n_replications=10,k=2,graph = FALSE)

## -----------------------------------------------------------------------------
summary(res_SIR_boot)

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res_SIR_boot,choice="size")

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res_SIR_boot,choice="selec_var")

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res_SIR_boot,choice="lambdas_replic")

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res_SIR_boot,choice="coefs_b")

## -----------------------------------------------------------------------------
wine <-  read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data", header = FALSE) 
colnames(wine) <- c('Type', 'Alcohol', 'Malic', 'Ash', 'Alcalinity', 'Magnesium', 
                    'Phenols', 'Flavanoids', 'Nonflavanoids', 'Proanthocyanins', 
                    'Color', 'Hue', 'Dilution', 'Proline')

# Extract the response variable
Y <- wine$Nonflavanoids
# Remove the response variable (Nonflavanoids) and the class information (type of cultivars)
X <- wine[, -which(names(wine) %in% c("Type","Nonflavanoids"))] 

head(cbind(Y,X),3)
print(dim(X))

## -----------------------------------------------------------------------------
X = scale(X)
Y = scale(Y)

## -----------------------------------------------------------------------------
res1 = SIR_threshold_opt(Y=Y, X=X, H=5, n_lambda=300, thresholding="soft", graph=FALSE)
summary(res1)

## ----fig.align='center',fig.width=5,fig.asp=1,out.width="300px"---------------
plot(res1,choice="estim_ind")

## ----fig.asp=1,fig.show='hold',fig.width=8,out.width="400px",fig.align='center',fig.cap=' '----
plot(res1,choice="cos2_selec")
plot(res1,choice="regul_path")

## -----------------------------------------------------------------------------
# To lighten the build of the vignette, the result of this command was saved into a RData file.
#res2 = SIR_threshold_bootstrap(Y=Y, X=X, H=5, n_lambda=200, thresholding="hard", n_replications = 500 , graph=FALSE)
load("../R/sysdata.rda") # load res2
summary(res2)

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res2,choice="size")

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res2,choice="selec_var")

## ----fig.align='center',fig.width=9,fig.asp=1,out.width="400px"---------------
plot(res2,choice="coefs_b")

## -----------------------------------------------------------------------------
res3 = SIR_threshold(Y=Y, X=X, H=5, lambda = 0.071, thresholding="hard", graph=FALSE)
summary(res3)

## -----------------------------------------------------------------------------
res4 = SIR(Y=res3$Y, X=res3$X_reduced, H=res3$H, graph = FALSE)
summary(res4)

## ----fig.align='center',fig.width=5,fig.asp=1,out.width="300px"---------------
plot(res4,choice="estim_ind")

## -----------------------------------------------------------------------------
summary(lm(Y~X))

## -----------------------------------------------------------------------------
summary(lm(Y~Flavanoids+Ash,data=wine))$r.squared

## -----------------------------------------------------------------------------
summary(lm(Y~Ash+Magnesium+Phenols+Flavanoids+Hue+Dilution,data=wine))$r.squared

