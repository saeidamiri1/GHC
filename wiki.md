# General Hybrid Clustering
Clustering is an unsupervised technique to find underlying structure in a dataset by grouping
data points into subsets that are as homogeneous as possible, clustering is a widely used unsupervised technique for identifying natural classes within a set of data. [Amiri et al. (2018)](https://github.com/saeidamiri1/GHC/blob/master/GHCsource/manuscript/manuscript-vr3.pdf) proposed a clustering technique for general clustering problems including those that have non-convex clusters. The proposed is fully nonparametric and it generates clusters for a given desired number of clusters K. They also discussed estimating the size of cluster.


## Contents
1. [Techniques](#techniques)
   - [Clustering](#clustering)
   - [Size of Cluster](#size-of-cluster)
2. [Data set](#datasets)
   - [Spiral data](#spiral-data)
   - [GARBER data](#garber-data)
   - [Wheat metabolomics data](#wheat-metabolomics-data)
   - [Scales data](#scales-data)
   - [Flame data](#flame-data)
3. [How to run the proposed hybrid clustering](how-to-run)
   - [Upload Source](#upload-source)
   - [Prepare data set](#prepare-data-set)
   - [Run Cluster](#run-cluster)
   - [Size of clusters](#size-of-clusters)
4. [GARBER data](#garber-data)
    - [kmeans](#kmeans)
    - [SHC](#shc)
    - [Spectral clustering](#spectral-clustering)
    - [CT](#ct)
    - [Chameleon (CHA)](#chameleon-cha)
    - [Trimmed clustering](#trimmed-clustering)
    - [Fuzzy clustering](#fuzzy-clustering)
    - [tSNE](#tsne)

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

The source is also available on [Github](https://github.com/saeidamiri1/GHC/tree/master/dataset/metabolicwheat).

```
whme<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/metabolicwheat/datasheet1.csv",sep='',header = TRUE)
```


#### Scales data
We generated many simulated convex data sets with  very different scales in one dimension (vertical) but similar scales on another (horizontal). One of them is Github and accessible via the following scripts   

```
scale0<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/scales0.csv",sep=",",header=TRUE)
scaled<-as.matrix(scale0[,2:3])
scalel<-scale0[,1]
```

#### FLAME data
Fu and Medico (2007) developed a fuzzy clustering technique for DNA microarray data which they considered on the test data FLAME, the data is accessible via the following script.

```
flame0<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/flame0.csv",sep=",",header=TRUE)
flamed<-as.matrix(flame0[,2:3])
flamel<-flame0[,4]
```

# How to run
### Upload Source
The source of codes are available in GitHub and using the following code can be uploaded in R
```
source('https://raw.githubusercontent.com/saeidamiri1/GHC/master/codes/SHC.R')
```

Also load the following libraries which run the computations in parallel,

```
library("foreach")
library("doParallel")
```

###  Prepare data set
To describe the codes, we used the spiral data,  the following script load source, dataset and plot it.  Use your data set with **matrix format**.

```
> spiral0<-read.csv("https://raw.githubusercontent.com/saeidamiri1/GHC/master/dataset/spiral0.csv",sep=",",header=TRUE)
> spiral<-as.matrix(spiral0[,2:3])
> head(spiral)
        V1   V2
[1,] 31.95 7.95
[2,] 31.15 7.30
[3,] 30.45 6.65
[4,] 29.70 6.00
[5,] 28.90 5.55
[6,] 28.05 5.00
```

### Run Cluster
Once the data and the codes are loaded in R, the clustering can be obtained using the following script

```
> knmin0<-2
> knmax0<-floor(dim(spiral)[1]/5)
> knmax0
[1] 62
> CLUS<-SHC(as.matrix(spiral),3,B=200,knmin=knmin0,knmax=knmax0)
```

The dendrogram can be also plotted,
```
> # plot the dendrogram
> plot(hclust(CLUS[[1]],method="single"),h=-1)
```

<img src="https://github.com/saeidamiri1/GHC/blob/master/images/Rplot01.jpeg" width="500">


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
plot(spiral,col=CLUS[[2]])
```
<img src="https://github.com/saeidamiri1/GHC/blob/master/images/Rplot02.jpeg" width="300">


It is of interest to run the proposed method for a cluster size of #4,  the following shows the codes and the clusters,

```
> CLUS<-SHC(Spiral,4,B=200,knmin=knmin0,knmax=knmax0)
> plot(spiral,col=CLUS[[2]])
```

<img src="https://github.com/saeidamiri1/GHC/blob/master/images/Rplot03.jpeg" width="300">


###  Size of clusters
 The following script shows the function of estimating the size of clusters,

```
> KCLUS<-EK(spiral,B=200,knmin=knmin0,knmax=knmax0)
```

```
# plot the dendrogram
> plot(hclust(KCLUS[[1]],method="single"),h=-1)
```
<img src="https://github.com/saeidamiri1/GHC/blob/master/images/Rplot04.jpeg" width="500">


```
> # print the assigned clusters to observation
> print(KCLUS[[2]])
[1] 3
```

##  GARBER data
To explain different  clustering methods, we run the clustering with [GARBER data](#garber-data).

```
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

#### kmeans

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


#### Spectral clustering
We also used spectral clustering, see von Luxburg et al. (2008)

```
library("kernlab")
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
> AI(cutree(hybridHclust(garber),6),cagarber)
[1] 0.6338028
```

#### Chameleon (CHA)
The CHAMELEON (CHA) clustering, Karypis et al. (1999),  is implemented in a software entitled ["cluto"](http://glaros.dtc.umn.edu/gkhome/cluto/cluto/overview), we used cluto-2.1.2 for the computation.  To run the clustering download the software from http://glaros.dtc.umn.edu/gkhome/cluto/cluto/download. The manual is inside the software, we put is in the Github as well [Chameleon manual]((https://github.com/saeidamiri1/GHC/blob/master/codes/clustomanual.pdf).



#### Trimmed clustering
garcia et al. (2008) consider a robustified form of clustering called trimmed clustering (TC), the central idea is that
the true clustering corresponds to a collection of normal distributions contaminated by outliers.

```
>library("tclust")
> AI(tkmeans(garber, k = 6, alpha=0.0)$cluster,cagarber)
[1] 0.6901408

```

#### Fuzzy clustering
To  run the fuzzy clustering,  we use, FCLUST.

```
> library("fclust")
>Cf<-AI(FKM(garber,k=6)$clus[,1],cagarber)
>Cf
[1] 0.3123239
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


## References
Amiri, S., Clarke, B, Clarke, J. & Koepke, H.A. (2018). A General Hybrid Clustering Technique. Accepted in Journal of Computational and Graphical Statistics. ([pdf](https://github.com/saeidamiri1/GHC/blob/master/manuscript/manuscript-vr3.pdf), [journal](https://www.tandfonline.com/toc/ucgs20/current))

Chipman, H. and R. Tibshirani (2006). Hybrid hierarchical clustering with applications to microarray data. Biostatistics 7(2), 286–301.

Karypis, G., E.-H. Han, and V. Kumar (1999). Chameleon: Hierarchical clustering using dynamic modeling. Computer 32 (8), 68–75. http://glaros.dtc.umn.edu/gkhome/ cluto/cluto/overview Accessed: 2018-07-20.

Fu, L. and E. Medico (2007). FLAME, a novel fuzzy clustering method for the analysis of dna microarray data. BMC Bioinformatics 8, 3.


Garber, M., O. Troyanskaya, K. Schluens, S. Petersen, Z. Thaesler, M. Pacyna-Gengelbach, Van De Rijn, G. Rosen, C. Perou, R. Whyte, Alman, D. Brown P, Botstein, and I. Petersen (2001). Diversity of gene expression in adenocarcinoma of the lung. Proceedings of the National Academy of Sciences 98(24), 13784–13789.

García-Escudero, L. A., Gordaliza, A., Matrán, C., & Mayo-Iscar, A. (2008). A general trimming approach to robust cluster analysis. The Annals of Statistics, 36(3), 1324-1345.


Maaten, L. V. D., & Hinton, G. (2008). Visualizing data using t-SNE. Journal of machine learning research, 9(Nov), 2579-2605.

Kessler, N., Bonte, A., Albaum, S. P., Mäder, P., Messmer, M., Goesmann, A., ... & Nattkemper, T. W. (2015). Learning to classify organic and conventional wheat–a machine learning driven approach using the MeltDB 2.0 metabolomics analysis platform. Frontiers in bioengineering and biotechnology, 3, 35.

von Luxburg, U., M. Belkin, and O. Bousquet (2008). Consistency of spectral clustering. Ann. Stat. 36, 555–586.

**[⬆ back to top](#contents)**
