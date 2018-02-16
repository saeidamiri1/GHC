source('https://raw.githubusercontent.com/saeidamiri1/GHC/master/SHC.R')
library("foreach")
library("doParallel")


datasource <- "https://github.com/saeidamiri1/GHC/blob/master/SPIRAL.RData?raw=true"
 load(url(datasource))
 Spiral
 plot(Spiral)

 knmin0<-2
 knmax0<-floor(dim(Spiral)[1]/5)
 knmax0
 CLUS<-SHC(Spiral,3,B=200,knmin=knmin0,knmax=knmax0)
 
# plot the dendrogram 
plot(hclust(CLUS[[1]],method="single"),h=-1)

# print the assigned clusters to observation
print(CLUS[[2]])

# plot the data with the assigned clusters
plot(Spiral,col=CLUS[[2]])
 



KCLUS<-EK(Spiral,B=200,knmin=knmin0,knmax=knmax0)
# plot the dendrogram 
plot(hclust(KCLUS[[1]],method="single"),h=-1)

# print the assigned clusters to observation
print(KCLUS[[2]])
