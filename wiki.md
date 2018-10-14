# General Hybrid Clustering
Clustering is an unsupervised technique to find underlying structure in a dataset by grouping
data points into subsets that are as homogeneous as possible, clustering is a widely used unsupervised technique for identifying natural classes within a set of data. [Amiri et al. (2018)](https://github.com/saeidamiri1/GHC/blob/master/manuscript/manuscript-vr3.pdf) proposed a clustering technique for general clustering problems including those that have non-convex clusters. The proposed is fully nonparametric and it generates clusters for a given desired number of clusters K. They also discussed estimating the size of cluster.


## Contents
1. [Techniques](#techniques)
   - [Clustering](#clustering)
   - [Size of Cluster](#size-of-cluster)
2. [Data set](#datasets)
   - [spiral data](#)
   - [GARBER data](#)
   - [wheat metabolomics data](#)
   - [Scales](#)
   - [FLAME data](#)

3. [Working with data]
   - [Source](#source)
   - [Data set](#data-set)
   - [Run Cluster](#run-cluster)
   - [Size of clusters](#size-of-clusters)

4. [GARBER data](#GARBER data)
5. [References](#references)

# Techniques
 We developed algorithms for clustering and estimating the size of cluster which are explained in the belows.
## Clustering
The proposed clustering method is referred to as Stabilized Hybrid Clustering (SHC) and its steps is presented in Algorithm 1,

<img src="https://github.com/saeidamiri1/GHC/blob/master/images/algorithm1.png" width="800">

Algorithm 1 is implemented in R,

```
SHC(x,K,B=200,knmin,knmax)
```

The arguments are: ```x``` is the observation, use the R's matrix format. ```B``` is number of run to get a stabilized clusters, we used B=200 in our computations. Concerning ```B```, run the code with different Bs and if you see huge different in result, increase the number of iterations. knmin and knmax are the minimum and maximum size of cluster to get the stabilized clustering. We used knmin=2 and knmax=n/5, where n is the sample size.
This ```SHC()``` provides the distance matrix and the predicted cluster.


## Size of Cluster
[Amiri et al. (2018)](https://github.com/saeidamiri1/GHC/blob/master/manuscript/manuscript-vr3.pdf) also discussed a technique to estimate the size of clusters, it is presented in Algorithm 2,

<img src="https://github.com/saeidamiri1/GHC/blob/master/images/algorithm2.png" width="800">

Algorithm 2 is implemented in R,
```
EK(observation,B=200,knmin,knmax)
```


# Data sets
#### Spiral data
The spiral data which is a non-convex data, the data can be upload via the following script, we used this data to explain the proposed algorithm.
```
spiral0<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/spiral0.csv",sep=",",header=TRUE)
spiral<-spiral0[,2:3]
plot(spiral)
```
<img src="https://github.com/saeidamiri1/GHC/blob/master/images/spiral.jpeg" width="300">



####  GARBER data
To study the performance of the SHC’s with high dimensional data, we used the microarray data from Garber et al. (2001). The data are the 916-dimensional gene expression profiles for lung tissue from n = 72 subjects. Of these, five subjects were normal and 67 had lung tumors. The classification of the tumors into 6 classes (plus normal) was done by a pathologist giving seven classes total.

```
library("pvclust")
data(lung)
attach(lung)

 garber<-t(lung)
 garber<-garber[-c(1,20),]

 for(i in 1:(dim(garber)[2])){
   garber[is.na(garber[,i]),i]<-mean(garber[,i],na.rm=T)
 }
 garber<-as.matrix(garber)

 # extract the tru label
 row1<-grep("Adeno",row.names(garber), perl=TRUE, value=FALSE)
 row2<-grep("normal",row.names(garber), perl=TRUE, value=FALSE)
 row3<-grep("SCLC",row.names(garber), perl=TRUE, value=FALSE)
 row4<-grep("SCC",row.names(garber), perl=TRUE, value=FALSE)
 row5<-grep("node",row.names(garber), perl=TRUE, value=FALSE)
 row6<-grep("LCLC",row.names(garber), perl=TRUE, value=FALSE)

 cagarber<-NULL
 cagarber[row1]<-1
 cagarber[row2]<-2
 cagarber[row3]<-3
 cagarber[row4]<-4
 cagarber[row5]<-5
 cagarber[row6]<-6

```

#### Wheat metabolomics data
The other dataset that we used is about Wheat metabolomics data, see Kessler et al. (2015). The dataset has 313 observations with 37 quantitative variables, and we consider ”Variety” (Antonius, CCP, Caphorn, DJ, MC2,Probus ,RdB ,R, Sandomir, Scaro, Titlis) as true label to evaluate the clustering techniques.  This data is accesible via the following script.


```
whme<-read.csv("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4371749/bin/Data_Sheet_1.CSV",sep='',header = TRUE)
```

The source is also availabe on [github](https://github.com/saeidamiri1/GHC/tree/master/dataset/metabolicwheat).


#### Scales
read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/flame0.csv",sep=",")

```
scale0<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/scale0.csv",sep=",",header=TRUE)
scaled<-scale0[,2:3]
scalel<-scale0[,1]
```

#### FLAME data
Fu and Medico (2007) developed a fuzzy clustering technique for DNA microarray data which they considered on the test data FLAME, the data is accessible via the following script.

```
flame0<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/flame0.csv",sep=",",header=TRUE)
flamed<-flame0[,2:3]
flamel<-flame0[,4]
```




#  Manual
### Source
The source of codes are available in GitHub and using the following code can be uploaded in R
```
> source('https://raw.githubusercontent.com/saeidamiri1/GHC/codes/master/SHC.R')
```

Also load the following libraries which run the computations in parallel,

```
> library("foreach")
> library("doParallel")
```

###  Prepare data set
To describe the codes, we used the spiral data,  the following script load source, dataset and plot it,

```
> source('https://raw.githubusercontent.com/saeidamiri1/GHC/codes/master/SHC.R')
> library("foreach")
> library("doParallel")

> datasource <- "https://github.com/saeidamiri1/GHC/datasets/blob/master/SPIRAL.RData?raw=true"
> load(url(datasource))
> Spiral
> plot(Spiral)
```


### Run Cluster
Once the data and the codes are loaded in R, the clustering can be obtained using the following script

```
> knmin0<-2
> knmax0<-floor(dim(Spiral)[1]/5)
> knmax0
[1] 62
> CLUS<-SHC(Spiral,3,B=200,knmin=knmin0,knmax=knmax0)
```

The dendrogram can be also plotted,
```
> # plot the dendrogram
> plot(hclust(CLUS[[1]],method="single"),h=-1)
```

<img src="https://github.com/saeidamiri1/GHC/image/blob/master/Rplot01.jpeg" width="500">


The predicted clusters are also available,

```
> # print the assigned clusters to observation
> print(CLUS[[2]])
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [42] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [83] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
[124] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
[165] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
[206] 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
[247] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
[288] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
```

```
# plot the data with the assigned clusters
plot(Spiral,col=CLUS[[2]])
```
<img src="https://github.com/saeidamiri1/GHC/image/blob/master/Rplot02.jpeg" width="300">


It is of interest to run the proposed method for a cluster of #4,  the following shows the codes and the clusters,

```
> CLUS<-SHC(Spiral,4,B=200,knmin=knmin0,knmax=knmax0)
> plot(Spiral,col=CLUS[[2]])
```

<img src="https://github.com/saeidamiri1/GHC/image/blob/master/Rplot03.jpeg" width="300">



###  Size of clusters
 The following script shows the function of estimating the size of clusters,

```
> KCLUS<-EK(Spiral,B=200,knmin=knmin0,knmax=knmax0)
[1] 3
```

```
# plot the dendrogram
> plot(hclust(KCLUS[[1]],method="single"),h=-1)
```
<img src="https://github.com/saeidamiri1/GHC/image/blob/master/Rplot04.jpeg" width="500">


```
> # print the assigned clusters to observation
> print(KCLUS[[2]])
[1] 3
```

##  GARBER data
To study the performance of the SHC’s with high dimensional data, we used the microarray data from Garber et al. (2001). The data are the 916-dimensional gene expression profiles for lung tissue from n = 72 subjects. Of these, five subjects were normal and 67 had lung tumors. The classification of the tumors into 6 classes (plus normal) was done by a pathologist giving seven classes total.

```
>
> ########################
> library("pvclust")
> data(lung)
> attach(lung)
....

> garber<-t(lung)
> garber<-garber[-c(1,20),]
>
> for(i in 1:(dim(garber)[2])){
+   garber[is.na(garber[,i]),i]<-mean(garber[,i],na.rm=T)
+ }
> garber<-as.matrix(garber)
>
> # extract the tru label
> row1<-grep("Adeno",row.names(garber), perl=TRUE, value=FALSE)
> row2<-grep("normal",row.names(garber), perl=TRUE, value=FALSE)
> row3<-grep("SCLC",row.names(garber), perl=TRUE, value=FALSE)
> row4<-grep("SCC",row.names(garber), perl=TRUE, value=FALSE)
> row5<-grep("node",row.names(garber), perl=TRUE, value=FALSE)
> row6<-grep("LCLC",row.names(garber), perl=TRUE, value=FALSE)
>
>
> cagarber<-NULL
> cagarber[row1]<-1
> cagarber[row2]<-2
> cagarber[row3]<-3
> cagarber[row4]<-4
> cagarber[row5]<-5
> cagarber[row6]<-6
>
```
#### kmean
```
> Ckm<-NULL
> for (i in 1:200){
+   Ct0<-kmeans(garber,6)
+   Ckm[i]<-AI(Ct0$cluster,cagarber)  
+ }
> mean(Ckm);sd(Ckm)
[1] 0.6766901
[1] 0.05749914
```
#### SHC

```
>
> knmin0<-2
> knmax0<-floor(dim(garber)[1]/5)
> knmax0
[1] 14
>
>
> ClasPred<-SHC(garber,6,B=200,knmin=knmin0,knmax=knmax0)
>
> ## accuracy index (AI)
> library('mclust')
> AI<-function(trueclu,obsclus)  1-classError(trueclu, obsclus)$errorRate
>
> AI(ClasPred[[2]],cagarber)
[1] 0.8169014
```

#### tSNE

Maaten and Hinton (2008) considered the t-Distributed Stochastic Neighbor Embedding (t-SNE) that is a technique for dimensionality reduction which was developed to for the visualization of high-dimensional datasets, the following shows the AI of using  t-SNE with SHC,

```
>library(tsne)

> Ctsh<-Ctkm<-NULL
> for (i in 1:40){
+   tsne_garber = tsne(garber, perplexity=50)
+   Ct1<-kmeans(tsne_garber,6)
+   Ct2<-SHC(tsne_garber,6,B=200,knmin=knmin0,knmax=knmax0)
+   Ctkm[i]<-AI(Ct1$cluster,cagarber)  
+   Ctsh[i]<-AI(Ct2[[2]],cagarber)
+ }
......

> mean(Ctkm);sd(Ctkm)
[1] 0.3528169
[1] 0.03649867
>
> mean(Ctsh);sd(Ctsh)
[1] 0.5619718
[1] 0.02339459
```

#### Spectral

```
> Ctsp<-NULL
> for (i in 1:100){
+   Ctsp[i]<-avcorrec2(c(specc(garber, centers=6)),cagarber)
+ }
> mean(Ctsp);sd(Ctsp)
[1] 0.734507
[1] 0.06801622
```

#### CT
CT is Hybrid hierarchical clustering using mutual clusters developed in Chipman and Tibshirani (2006),  Code implementing CT is in the R package hybridHclust.

```
> library("hybridHclust")
> avcorrec2(cutree(hybridHclust(garber),6),cagarber)
[1] 0.6338028
```

#### Fuzzy clustering

```
> library("fclust")
>Cf<-AI(FKM(garber,k=6)$clus[,1],cagarber)
>Cf
[1] 0.3123239

```







### References
Amiri, S., Clarke, B, Clarke, J. & Koepke, H.A. (2018). A General Hybrid Clustering Technique. Accepted in Journal of Computational and Graphical Statistics. ([pdf](https://github.com/saeidamiri1/GHC/blob/master/manuscript/manuscript-vr3.pdf), [journal](https://www.tandfonline.com/toc/ucgs20/current))

Chipman, H. and R. Tibshirani (2006). Hybrid hierarchical clustering with applications to microarray data. Biostatistics 7(2), 286–301.

Fu, L. and E. Medico (2007). FLAME, a novel fuzzy clustering method for the analysis of dna microarray data. BMC Bioinformatics 8, 3.


Garber, M., O. Troyanskaya, K. Schluens, S. Petersen, Z. Thaesler, M. Pacyna-Gengelbach, Van De Rijn, G. Rosen, C. Perou, R. Whyte, Alman, D. Brown P, Botstein, and I. Petersen (2001). Diversity of gene expression in adenocarcinoma of the lung. Proceedings of the National Academy of Sciences 98(24), 13784–13789.

Maaten, L. V. D., & Hinton, G. (2008). Visualizing data using t-SNE. Journal of machine learning research, 9(Nov), 2579-2605.

Kessler, N., Bonte, A., Albaum, S. P., Mäder, P., Messmer, M., Goesmann, A., ... & Nattkemper, T. W. (2015). Learning to classify organic and conventional wheat–a machine learning driven approach using the MeltDB 2.0 metabolomics analysis platform. Frontiers in bioengineering and biotechnology, 3, 35.
