library(SIRthresholded)

a = read.csv("other/ozone.txt",sep = "")
Y = data.matrix(a[1],rownames.force = NA)
X = data.matrix(a[2:11],rownames.force = NA)

