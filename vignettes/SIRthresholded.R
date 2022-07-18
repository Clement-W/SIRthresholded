## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo=TRUE,eval=TRUE,comment="#>")

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

## ----fig.align='center',fig.width=4,fig.asp=1---------------------------------
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

## ----fig.asp=1,fig.show='hold'------------------------------------------------
plot(res_SIR,choice="estim_ind")
plot(res_SIR,choice="eigvals")

## -----------------------------------------------------------------------------
res_STSIR = SIR_threshold_opt(Y=Y, X=X, H=10, n_lambda=100, thresholding="soft", graph=FALSE)

## -----------------------------------------------------------------------------
summary(res_STSIR)

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_STSIR,choice="cos2_selec")

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_STSIR,choice="opt_lambda")

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_STSIR,choice="regul_path")

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_STSIR,choice="estim_ind")

## ----fig.asp=1,fig.show='hold',fig.width=6------------------------------------
res_HTSIR = SIR_threshold_opt(Y=Y, X=X, H=10, n_lambda=100, thresholding="hard", graph=FALSE)

plot(res_HTSIR,choice="cos2_selec")
plot(res_HTSIR,choice="regul_path")

## -----------------------------------------------------------------------------
res_SIR_pstar = SIR(Y=res_STSIR$Y, X=res_STSIR$X_reduced, H=res_STSIR$H)
summary(res_SIR_pstar)

## ----fig.asp=1,fig.width=4,fig.align="center"---------------------------------
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

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_SIR_boot,choice="size")

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_SIR_boot,choice="selec_var")

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_SIR_boot,choice="lambdas_replic")

## ----fig.asp=1,fig.width=6,fig.align='center'---------------------------------
plot(res_SIR_boot,choice="coefs_b")

