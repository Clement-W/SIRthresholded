library(SIRthresholded)

a = read.csv("other/ozone.txt",sep = "")
Y = data.matrix(a[1],rownames.force = NA)
X = data.matrix(a[2:11],rownames.force = NA)

a=SIR_threshold_bootstrap(Y,X,n_replications = 100,thresholding = "soft",n_lambda = 500)

b = read.csv("other/uscrime.txt",sep="\t")
Y = data.matrix(b[16],rownames.force = NA)
X = data.matrix(b[1:15],rownames.force = NA)
