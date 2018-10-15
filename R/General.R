
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


