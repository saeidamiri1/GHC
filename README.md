# GHC: 
This supplementary includes the codes of 
 
Amiri, S., Clarke, B, Clarke, J. & Koepke, H.A. (2017). A General Hybrid Clustering Technique.

Any feedback is really appreciated, please report bugs, typos or any comments by sending an email, saeid.amiri1 atsign gmail.com. 



# Load Spiral data from Github
datasource <- "https://github.com/saeidamiri1/GHC/blob/master/SPIRAL.RData?raw=true"
load(url(datasource))
Spiral
plot(Spiral)

# Load the source from Github
source('https://raw.githubusercontent.com/saeidamiri1/GHC/master/SHC.R')


library("foreach")
library("doParallel")

# To find the cluster of size of 3
CLUS<-SHC(Spiral,3,B=200)
CLUS
plot(Spiral,col=CLUS)


# To estimate the size of clusters
EK(Spiral,B=200)

