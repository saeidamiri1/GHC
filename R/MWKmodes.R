MWKModes<-function(data,center,a,T1,T2){
  
  center=sort(center)
  r<-dim(data)[1]
  cc<-dim(data)[2]
  k<-length(center)
  modes=matrix(0,nrow=k,ncol=cc)
  OLDFSS=0;
  clen<-NULL
  attr<-matrix(0,nrow=r,ncol=cc)# extra added
  for (i in 1:cc){
    column=sort(unique(data[,i]))
    clen[i]=length(column)
    for (j in 1:clen[i]){
      attr[j,i]=column[j]
    }
  }
  mclen<-max(clen)# extra added
  attr<-attr[1:mclen,]# extra added
  qus=0;
  ar=dim(attr)[1];ac=dim(attr)[2]
  attr_count=matrix(0,nrow=r,ncol=cc)#zeros(ar,ac);
  for (i in 1:r){#for i=1:r
    for (j in 1:cc){#for j=1:cc
      w=c(which(attr[,j]==data[i,j],arr.ind=T))#w=find(attr(:,j)==data(i,j));
      x=w[1]#x=w(1);
      attr_count[x,j]=attr_count[x,j]+1;  #attr_count(x,j)=attr_count(x,j)+1;
    }#end
  }#end
  
  
  for (i in 1:cc){#for i=1:cc
    for (j in 1:(clen[i])){#for j=1:clen(i)
      x=length(which(data[,i]==attr[j,i]))#x=length(find(data(:,i)==attr(j,i)));
      qus=qus+(r-x)/k#qus=qus+(r-x)/k;
    }#end
  }#end
  
  
  modeSum<-matrix(1,nrow=k,ncol=cc)#modeSum=ones(k,cc);
  for(i in 1:k){#for i=1:k
    modes[i,]=data[center[i],]  #modes(i,:)=data(center(i),:);
  }#end
  
  t=0;
  new_class=matrix(0,nrow=1,ncol=r)#new_class=zeros(1,r);
  old_class=matrix(0,nrow=1,ncol=r)#old_class=zeros(1,r);
  weight<-array(0,dim=c(2,k,cc))#weight=zeros(2,k,cc);
  weight[2,,]=1/cc#weight(2,:,:)=1/cc;
  weight[1,,]=1/cc#weight(1,:,:)=1/cc;
  
  while(1==1){#while 1==1
    
    t=t+1;
    new_mode=array(0,dim=c(k,ar,cc))#new_mode=zeros(k,ar,cc);
    classsum=matrix(0,nrow=1,ncol=k)#classsum=zeros(1,k);
    NEWFSS=0;
    for(i in 1:r){#for i=1:r
      distance<-matrix(0,nrow=1,ncol=k)#distance=zeros(1,k);
      for (j in 1:k){#for j=1:k
        
        for(h in 1:cc){#for h=1:cc
          if(data[i,h]!=modes[j,h])#if ~isequal(data(i,h),modes(j,h))
          {distance[j]=distance[j]+1+(weight[1,j,h])^a}
          else{
            distance[j]=distance[j]+(weight[2,j,h])^a#;%*modeSum(j,h)
          }
        }
      }
      
      MinValue=min(distance)
      MinRow=which.min(distance)
      old_class[i]=MinRow;
      NEWFSS=NEWFSS+MinValue;
      classsum[MinRow]=classsum[MinRow]+1;
      for (e in 1:cc){#for e=1:cc
        w=(which(attr[,e]==data[i,e]))[1]#w=find(attr(:,e)==data(i,e),1);
        # %x=w(1);
        new_mode[MinRow,w,e]=new_mode[MinRow,w,e]+1#new_mode(MinRow,w,e)=new_mode(MinRow,w,e)+1;
      }
    }
    
    
    for(i in 1:k){#for i=1:k
      for(h in 1:cc){#for h=1:cc
        NEWFSS=NEWFSS+T1*weight[1,i,h]^a+T2*weight[2,i,h]^a;
      }
    }
    new_sum=matrix(0,nrow=k,ncol=cc)#zeros(k,cc);
    modeSum=matrix(0,nrow=k,ncol=cc)#zeros(k,cc);
    for (i in 1:k){#for i=1:k       
      s_sum=0; 
      z_sum=0;
      for(j in 1:cc){ #for j=1:cc     
        if (classsum[i]!=0){
          #[Maxe,Maxr]=max(new_mode(i,:,j));
          Maxe<-max(new_mode[i,,j])[1]
          Maxr<-which.max(new_mode[i,,j])[1]  
          #if classsum(i)~=0
          modes[i,j]=attr[Maxr,j]
          modeSum[i,j]=(classsum[i]-Maxe)+T1
          new_sum[i,j]=Maxe+T2}
        else{          
          #%modeSum(i,:)=ones(1,cc);
          new_sum[i,j]=1;
        }
        s_sum=s_sum+(1/new_sum[i,j])^(1/(a-1));
        z_sum=z_sum+(1/modeSum[i,j])^(1/(a-1));
      }
      s_z=which(modeSum[i,(1:cc)]==0) #s_z=find(modeSum(i,[1:cc])==0);
      s_zl=length(s_z)
      if (s_zl!=0){
        weight[1,i,]=matrix(0,nrow=1,ncol=cc)#zeros(1,cc);
        for (h in 1:s_zl){
          weight[1,i,s_z[h]]=1/s_zl;  
        }}#end}
      else{
        for(j in 1:cc){# j=1:cc 
          weight[1,i,j]=(1/modeSum[i,j])^(1/(a-1))/z_sum;
        }
      }
      
      for (j in 1:cc){ 
        weight[2,i,j]=(1/new_sum[i,j])^(1/(a-1))/s_sum
        #%weight(1,i,j)=(1/modeSum(i,j))^(1/(a-1))/z_sum;
      }
      
    }
    
    
    NEWFSS=round(NEWFSS*100000)/100000;
    
    if (OLDFSS==NEWFSS)
      break
    else{
      OLDFSS=NEWFSS;
      new_class=old_class;
    }  
  }
  return(new_class)
}

