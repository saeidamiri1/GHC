############# Oct-12-2018
############# https://github.com/saeidamiri1/GHC
############# saeid.amiri1@gmail.com 
###################################
################################# Clustering    
SHC<-function (x,K,B=200, knmin=2,knmax=floor(dim(x)[1]/5)){
  x<-cbind(x,rep(0,nrow(x)))
  Len<-dim(x)
  clusterO<-Len[2]
  b<-1
  # it is 4
  # knmin<-2;knmax<-min(25,floor(dim(x)[1]/5)-2)
  kn<-sample(c(knmin,knmax),1)
  #dd<-Hub2MQ(x,kn)
  #RE<-dd[,Len[2]]
  
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
  
  
  
  cl <- makeCluster(detectCores()-1) # create a cluster with n-1 cores
  registerDoParallel(cl) # register the cluster
  ens = foreach(i = 1:200, 
                .combine = "rbind", .export=c("Hub2MQ","distancematrix0","reorderf","distwo")) %dopar% {
                  fit1 <- RHub2MQ(x,kn,knmin,knmax)
                  fit1
                }
  stopCluster(cl) 
  
  
  #while(b<B){
  #  kn<-sample(c(knmin:knmax),1)
  #  dd<-Hub2MQ(x,kn)
  #  RE<-rbind(RE,dd[,Len[2]])    
  #  b<-b+1
  #}
  
  REDIST<-as.dist(distancematrixH(ens))
  REDISTT<-as.matrix(REDIST)
  
  hc <- hclust(REDIST,method = "single")
  zhh<-mean(Xsub(hc,K),na.rm=T)
  kstar<-length(unique(cutree(hc,h=zhh)))
  
  cc<-cutree(hc,kstar)
  kn<-K
  xl2<-x
  ni<-Len[1]
  for(i in 1:ni){
    xl2[i,clusterO]<-cc[i]
  }
  
  
  alpha0<-.05
  while(alpha0>0){
    xcc<-NULL
    for(j in unique(cc))   xcc[j]<-length(xl2[xl2[,clusterO]==j,clusterO])
    mino0<- which(xcc/dim(x)[1]<alpha0)
    main0<-setdiff((cc),mino0)
    
    if(length(main0)>(kn)) break
    alpha0<-alpha0/2
  }
  
  
  i<-1
  cc0<-NULL
  for(j in main0){
    cc0[cc==j]<-i
    i<-i+1
  }
  for(j in mino0){
    cc0[cc==j]<-i
    i<-i+1
  }
  
  #cc02<-cc[1:length(main0)]
  kmi<-length(mino0)
  kma2<-kma<-length(main0)
  
  cc1<-cc0
  #cc2<-setdiff((cc),main0)
  while(kma2>kn){
    #xl<-xl[,-clusterO] 
    xz<-list(NULL)
    for(i in unique(cc1)){
      xz[[i]]<-which(cc1==i)
    }
    kcc<-unique(cc1)
    XXXX<-distancematrix0SE3(REDISTT,c(1:kma2),xz)
    cc1[cc1==XXXX[2]]<-XXXX[1]
    cc1<-reorderfSE(cc1)
    kma2<-kma2-1
  }
  
  main02<-1:kma2
  mino02<-(kma2+1):(kma2+kmi)
  
  
  xz<-list(NULL)
  for(i in unique(cc1)){
    xz[[i]]<-which(cc1==i)
  } 
  
  
  xl2<-x
  ni<-Len[1]
  for(i in 1:ni){
    xl2[i,clusterO]<-cc1[i]
  }
  
  
  if(!length(mino0)==0){ 
    while( !length(mino02)==0){
      ind<-md<-NULL
      i0<-1
      for(i1 in mino02 ){
        d1<-NULL
        for(i2 in main02){
          d1<-c(d1,distwoA(REDISTT,xz[[i2]],xz[[i1]]))
        }
        ind[i0]<-which.min(d1)
        md[i0]<-min(d1,na.rm=T)
        i0<-i0+1
      }
      
      xl2[xl2[,clusterO]==mino02[which.min(md)],clusterO]=main02[ind[which.min(md)]]
      
      cc1[cc1==mino02[which.min(md)]]<-main02[ind[which.min(md)]]
      xz<-list(NULL)
      for(i in unique(cc1)){
        xz[[i]]<-which(cc1==i)
      }
      
      mino02<-setdiff(mino02,mino02[which.min(md)])  
      #if(length(mino1)==0) break
    }} 
  
  return(list(REDIST,xl2[,clusterO]))
}

############################### SIZE of cluster
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


#######################
######################
######################
Hub2MQ<-function(x,kn){
  clusterO<-dim(x)[2]
  zz<-sample(c(4:6),1)
  (cl <- kmeans(x[,-clusterO], floor(dim(x)[1]/zz),  algorithm = "MacQueen",iter.max = 50, nstart = 1))  
  xl<-cbind(x[,-clusterO],cl$cluster)
  xlc<-distancematrix0(xl)
  # xlcT<-distancematrix0T(xl)
  xlc<-as.dist(xlc)  
  hc <- hclust(xlc,method = "single")
  xk<-NULL
  if( kn>length(cl$size)) kn<-length(cl$size)-1
  cc<-cutree(hc,kn)
  
  xl2<-cbind(x,NA)
  ni2<-sort(unique(xl[,clusterO]))
  for(i in ni2){
    xl2[xl[,clusterO]==i,clusterO+1]<-cc[i]
  }
  return(xl2[,-clusterO])
}




###distancematrixHT<-function(data){
###  Len<-dim(data)
###  ss<-Len[2]
###  dismat<- matrix(NA,ncol=ss,ss)
###  for(i in 1:ss){
###    for(j in 1:ss){
###      dismat[i,j]<-(distancem(data[,i],data[,j]))
###    }
###  }
###  t(dismat)
###}



distancematrixH<-function(data){
  Len<-dim(data)
  ss<-Len[2]
  dismat<- matrix(NA,ncol=ss,ss)
  for(i in 1:ss){
    if (i==ss) break
    for(j in ((i+1):ss)){
      dismat[i,j]<-(distancem(data[,i],data[,j]))
    }
  }
  t(dismat)
}


distancem<-function(a,b){
  distan1<-0
  distan1 <- sum(a!= b)
  return(distan1)
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


distwoA2<-function(xlcT0,dat1,dat2){
  dat11<-unlist(dat1)
  dat22<-unlist(dat2)
  di00<-NULL
  ff<-1
  for(ii1 in dat11){
    for(ii2 in dat22){
      di00[ff]<-xlcT0[ii1,ii2]
      ff<-ff+1
    }}
  min(di00)
}


distwoA<-function(xlcT0,dat1,dat2){
  dat11<-unlist(dat1)
  dat22<-unlist(dat2)
  #d1<-dim(data1)
  #if(is.null(d1)) {data1<-t(as.matrix(data1));d1<-dim(data1)}
  #d2<-dim(data2)
  #if(is.null(d2)) {data2<-t(as.matrix(data2));d2<-dim(data2)}
  di00<-NULL#matrix(NA,nrow=d1*d2,ncol=2)
  ff<-1
  for(ii1 in dat11){
    for(ii2 in dat22){
      di00[ff]<-xlcT0[ii1,ii2]
      ff<-ff+1
    }}
  min(di00)
}


distancematrix0SE3<-function(XLL,W1,XZZ){
  ss<-length(W1)
  dismat<- matrix(NA,ncol=ss,ss)#array(NA, dim=c(1,Len[2]-1,ss,ss))
  for(i in W1){
    if (i==W1[ss]) break
    for(j in (W1[(which(W1==i)+1):ss])){
      dismat[i,j]<-distwoA2(XLL,XZZ[[i]],XZZ[[j]])
    }
  }
  t(dismat)
  eee<-which(t(dismat)==min(t(dismat),na.rm=T), arr.ind = TRUE)
  sort(eee[sample(dim(eee)[1],1),])
  
}

reorderfSE<-function(data){
  Len<-length(data)
  nk<-as.integer(names(table(data)))
  h<-0
  data0<-NULL
  for(i in nk){
    ij<- data==i
    h<-h+1
    data0[ij]<-h
  }
  return(data0)
}


CreatXCC<-function(xx){
  xxc<-list()
  len<-dim(xx)
  xxc[[1]]<-1
  for(l1 in 2:(len[1])){
    if(xx[l1,1]<0&xx[l1,2]<0) xxc[[l1]]<-l1
    if(xx[l1,1]<0&xx[l1,2]>0)  xxc[[l1]]<-c(l1,unlist(xxc[[xx[l1,2]]])) 
    if(xx[l1,1]>0&xx[l1,2]>0) xxc[[l1]]<-c(l1,unlist(xxc[[xx[l1,1]]]),unlist(xxc[[xx[l1,2]]]))
  }  
  xxc
}


Xsub<-function(hclust0,K0){
  xcc<-list()
  xx<-hclust0$merge
  len<-dim(xx)
  zh<-NULL
  xcc<-CreatXCC(xx)
  xc<-NULL
  heigh1<-hclust0$heigh
  inverse = function (f, lower = -100, upper = 100) {
    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  }
  square_inverse = inverse(function (x) length(unique((cutree(hclust0,h=x)))), min(heigh1), max(heigh1))
  Km<-which.min((heigh1-square_inverse(K0)$root)^2)
  for(l1 in (Km-1):(1)){
    if(l1 %in% xc) next
    if(xx[l1,1]<0&xx[l1,2]<0) {zh[l1]<-heigh1[l1]}
    if(xx[l1,1]<0&xx[l1,2]>0) {zh[l1]<-heigh1[l1];xc<-c(xc,xcc[[xx[l1,2]]]) }    
    if(xx[l1,1]>0&xx[l1,2]>0) {zh[l1]<-heigh1[l1];xc<-c(xc,xcc[[xx[l1,2]]],xcc[[xx[l1,1]]]) }
  }  
  return(zh)
}


