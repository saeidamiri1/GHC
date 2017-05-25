# Load Spiral data from Github
datasource <- "https://github.com/saeidamiri1/GHC/blob/master/SPIRAL.RData?raw=true"
load(url(datasource))
Spiral
plot(Spiral)

# Load the source from Github
source('https://raw.githubusercontent.com/saeidamiri1/GHC/master/SHC.R')


library("foreach")
library("doParallel")

# Find the cluster
CLUS<-SHC(Spiral,3,B=200)
CLUS
plot(Spiral,col=CLUS)



# Estimate the number of cluster
EK(Spiral,B=10)

