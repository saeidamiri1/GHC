############# Oct-13-2018
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
