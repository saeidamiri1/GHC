# Here the hidden functions are given.  
#### Hamming distance
distancem<-function(a,b){
n<-length(a)
distan1<-0
if(length(a)!=length(b)) distance = -1;
    for (i in 1:length(a)){
        distan1 <- distan1 + (a[i]!= b[i])
     }
return(distan1)
}


#distancem<-function(a,b){
#  n<-length(a)
#  distan1<-0
#  if(length(a)!=length(b)) {distance1 = -1}
#  else{distan1 <-sum(a-b!=0,na.rm = T)}
#  return(distan1)
#}


## Find Mode
Mode <- function(x) {
x<-sort(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
## Find Mode of column
Modv<-function(x){
Len<-dim(x)
re<-NULL
for (j in 1:Len[2]){
re[j]<-Mode(x[,j])
}
re}


#distancematrixH<-function(data){
  #please use with more than two dimention
  #data<-reorderf(data)
#  Len<-dim(data)
#  ss<-Len[2]
#  dismat<- matrix(NA,ncol=ss,ss)#array(NA, dim=c(1,Len[2]-1,ss,ss))
#  for(i in 1:ss){
#    if (i==ss) break
#    for(j in ((i+1):ss)){
#      dismat[i,j]<-(distancem(data[,i],data[,j]))
 #   }
#  }
#  t(dismat)
#}

######################################
# Original one
#HCSL<-function (x,kn0){
#	x<-x[,-(dim(x)[2])]
#REDIST<-as.dist(distancematrixH(t(x)))
#hc <- hclust(REDIST,method = "single")
#cutree(hc,k=kn0)
#}


#HCAL<-function (x,kn0){
#	x<-x[,-(dim(x)[2])]
#REDIST<-as.dist(distancematrixH(t(x)))
#hc <- hclust(REDIST,method = "average")
#cutree(hc,k=kn0)
#}


#HCCL<-function (x,kn0){
#	x<-x[,-(dim(x)[2])]
#REDIST<-as.dist(distancematrixH(t(x)))
#hc <- hclust(REDIST,method = "complete")
#cutree(hc,k=kn0)
#}

##########################
#ensemblecluobse<-function(x,B){
#  Len<-dim(x)
#  clusterO<-Len[2]
#  b<-1
#  knm<-ifelse(Len[1]>100, 70, Len[1]/2)
#  kn<- sample(c(2,knm),1)
#  dd<-HCAL(x,kn)
#  RE<-dd
#  while(b<B){
#    kn<-sample(c(2:knm),1)
#    dd<-HCAL(x,kn)
#    RE<-rbind(RE,dd)    
#    b<-b+1
#  }
#  return(RE)
#}
###############################################
######################
##### For high dimensional
###############################################
hammingD<-DistM<-function(Xdat){
  Len<-dim(Xdat)
  ss<-Len[1]
  dismat<- matrix(1,ncol=ss,ss)#array(NA, dim=c(1,Len[2]-1,ss,ss))
  for(i in 1:ss){
    if (i==ss) break
    for(j in (i:ss)){
      dismat[i,j]<-mean(Xdat[i,]==Xdat[j,],na.rm=T)#(distancem(data[,i],data[,j]))
    }
  }
  return(1-t(dismat))
}

#hammingD<-function(Xdat){
#  Len<-dim(Xdat)
#  ss<-Len[1]
#  dismat<- matrix(1,ncol=ss,ss)#array(NA, dim=c(1,Len[2]-1,ss,ss))
#  for(i in 1:ss){
#    if (i==ss) break
#    for(j in (i:ss)){
#      dismat[i,j]<-mean(Xdat[i,]==Xdat[j,],na.rm=T)#(distancem(data[,i],data[,j]))
#    }
#  }
#  return(1-t(dismat))
#}

###############################################
HCALB<-function (x,kn0){
  #x<-x[,-(dim(x)[2])]
  REDIST<-as.dist(DistM((x)))
  hc <- hclust(REDIST,method = "average")
  cutree(hc,k=kn0)
}

