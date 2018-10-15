############# Oct-12-2018
############# https://github.com/saeidamiri1/GHC
############# saeid.amiri1@gmail.com
###################################
############################### SIZE of cluster
EK<-function (x,B=200, knmin=2,knmax=floor(dim(x)[1]/5)){
  x<-cbind(x,rep(0,nrow(x)))
  Len<-dim(x)
  clusterO<-Len[2]
  b<-1
  #knmin<-2;knmax<-min(25,floor(dim(x)[1]/5)-2)
  kn<-sample(c(knmin,knmax),1)
  dd<-Hub2MQ(x,kn)
  RE<-dd[,Len[2]]

  #####
  ####
  RHub2MQ<-function(x,kn,knmin,knmax){
    kn<-sample(c(knmin:knmax),1)
    dd<-Hub2MQ(x,kn)
    Len<-dim(x)
    return(dd[,Len[2]])
  }

  distancematrix0<-function(data){
    data<-reorderf(data)
    Len<-dim(data)
    nk<-as.integer(names(table(data[,Len[2]])))
    ss<-length(unique(data[,Len[2]]))
    dismat<- matrix(NA,ncol=ss,ss)#array(NA, dim=c(1,Len[2]-1,ss,ss))
    for(i in 1:ss){
      if (i==ss) break
      for(j in ((i+1):ss)){
        ij0<-data[,Len[2]]==i
        ij1<-data[,Len[2]]==j
        dismat[i,j]<-(distwo(data[ij0,-Len[2]],data[ij1,-Len[2]]))
      }
    }
    t(dismat)
  }

  reorderf<-function(data){
    Len<-dim(data)
    nk<-as.integer(names(table(data[,Len[2]])))
    h<-0
    data0<-cbind(data,NA)
    for(i in nk){
      ij<-data[,Len[2]]==i
      h<-h+1
      data0[ij,Len[2]+1]<-h
    }
    return(data0[,-Len[2]])
  }

  distwo<-function(data1,data2){
    d1<-dim(data1)
    if(is.null(d1)) {data1<-t(as.matrix(data1));d1<-dim(data1)}
    d2<-dim(data2)
    if(is.null(d2)) {data2<-t(as.matrix(data2));d2<-dim(data2)}
    di0<-NULL
    ff<-1
    for(i in 1:d1[1]){
      for(j in 1:d2[1]){
        di0[ff]<-mean((data1[i,]-data2[j,])^2)
        ff<-ff+1
      }
    }
    quantile(di0,probs=.2)
  }



  cl <- makeCluster(detectCores()-1) # create a cluster with 2 cores
  registerDoParallel(cl) # register the cluster
  ens = foreach(i = 1:B,
                .combine = "rbind", .export=c("Hub2MQ","distancematrix0","reorderf","distwo")) %dopar% {
                  fit1 <- RHub2MQ(x,kn,knmin,knmax)
                  fit1
                }
  stopCluster(cl)


  REDIST<-as.dist(distancematrixH(ens))
  hclustM <- hclust(REDIST, method = "single")
  cutValue <- hclustM$height[which.max(diff(hclustM$height))]
  ee<-(cutree(hclustM, h = cutValue))
  ee0<-length(unique(ee))
  #ee<-(cutree(hclustM, h = cutValue))
  idn<-as.numeric(names(table(ee)))[table(ee)/length(ee)<.009]
  eeNA<-NULL
  for(i in idn){
    eeNA<-c(eeNA,which(ee==i))
  }


  if(length(eeNA)!=0){
    SEMAX<-sort(which(diff(hclustM$height)==sort(diff(hclustM$height), decreasing = TRUE)[2]),decreasing = TRUE)[1]
    cutValue <- hclustM$height[SEMAX]
    ee2<-(cutree(hclustM, h = cutValue))
    ee0<-mean(c(sum(table(ee2)/length(ee2)>.009),length(unique(ee))))
  }

  return(list(REDIST,ee0))

}


