data<-read.csv("new.csv",fileEncoding="latin1")
library(readr)
data$floor<-parse_number(data$floor)
data<-data[!is.na(data$floor),]
data<-data[!is.na(data$buildingType),]

data$tradeTime<-as.Date(data$tradeTime,format = "%m/%d/%y")

#save(data,file = "full beijing.RData")
#data<-load("full beijing.RData")

data1<-data[data$tradeTime>="2016-08-01"&data$tradeTime<="2016-09-30", ]
sum(!duplicated(data1[,1:2]))
data1<-data1[data1$totalPrice>=100,]
sum(!duplicated(data1[,1:2]))

data<-data.frame(index=1:dim(data1)[1],data1)
data<-data[order(data$Lng),]
coord<-data[,c("Lng","Lat")]

data$bathRoom<-as.numeric(data$bathRoom)
data$drawingRoom<-as.numeric(data$drawingRoom)
data$renovationCondition<-as.numeric(data$renovationCondition)
data$livingRoom<-as.numeric(data$livingRoom)


require(dplyr)
require(sp)
unique.coord<-rename(count(coord, coord[,1], coord[,2]), Freq = n)
unique.coord<-unique.coord[order(unique.coord[,1]),]
unique.coord<-data.frame(unique.coord, location=10000+(1:dim(unique.coord)[1]))
data$location<-rep(unique.coord$location,unique.coord$Freq)
colnames(unique.coord)[1:2]<-c("Lng","Lat")

loc.train<-sample(unique.coord$location,size = dim(unique.coord)[1]*0.8)
test.coord<-subset(unique.coord, !location%in%loc.train)
train.coord<-subset(unique.coord, location%in%loc.train)
train<-subset(data,location%in%loc.train)
test<-subset(data,!location%in%loc.train)


#save(test,test.coord,train,train.coord,loc.train,unique.coord,file = "beijing.RData")
load("beijing.RData")
require(sp)

D<-spDists(as.matrix(unique.coord[,1:2]),longlat = T)
quantile(D[lower.tri(D)],c(0.01,0.25,0.5,0.75,0.99))
summary(D[lower.tri(D)]) #U(0.03, 100)

M=20 #number of neighbors
NN_ind<-matrix(rep(0, (dim(train.coord)[1]-1)*M),nrow = dim(train.coord)[1]-1,ncol = M)
NN_distM<-matrix(rep(0, (dim(train.coord)[1]-1)*M*(M-1)/2),nrow = dim(train.coord)[1]-1,ncol = M*(M-1)/2)
NN_dist<-matrix(rep(0, (dim(train.coord)[1]-1)*M),nrow = dim(train.coord)[1]-1,ncol = M)

for (i in 2:dim(train.coord)[1]) {
  distance<-spDistsN1(as.matrix(train.coord[1:(i-1),1:2]),as.matrix(train.coord[i,1:2]),longlat = T)
  near.distance<-sort(distance,decreasing = FALSE,index.return=TRUE)
  dist.NN<-near.distance$x[1:min(i-1,M)] #dist to nearest neighbors
  ind.NN<-near.distance$ix[1:min(i-1,M)] #index of nearest neighbors
  distM.NN<-spDists(as.matrix(train.coord[ind.NN,1:2]),longlat = T)
  
  
  NN_ind[i-1,1:min(i-1,M)]<-ind.NN
  if(i>=3){
  NN_distM[i-1,1:(min(i-1,M)*(min(i-1,M)-1)/2)]<-distM.NN[lower.tri(distM.NN,diag = FALSE)]
  }
  NN_dist[i-1,1:min(i-1,M)]<-dist.NN
}

neighbor.train<-list(NN_ind=NN_ind, NN_distM=NN_distM, NN_dist=NN_dist)
require(ggplot2)
ggplot(data = unique.coord,aes(x=Lng,y=Lat))+geom_point()

######
###### Construct NNGP Prior within INLA
######

inla.rgeneric.nngp.model <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL) {
  
  envir = parent.env(environment())
  #Internal function
  
  interpret.theta <- function() {
    return(
      list(sigma2 = exp(-theta[1]),
           phis = (a+b*exp(theta[2]))/(1+exp(theta[2])))
    )
  }
  
  graph = function() {
    require(Matrix)
    
    return (Q())
  }
  
  Q = function() {
    require(Matrix)
    require(parallel)
    param = interpret.theta()
    
    sigma2<-param$sigma2
    phis<-param$phis
    
    NN_ind<-neighbor.train$NN_ind
    NN_ind<-rbind(rep(0,dim(NN_ind)[2]),NN_ind)
    B<-list()
    F<-list()
    for (j in 1:dim(NN_ind)[1]) {
      
      if(j==1){
        B[[j]]<-0
        F[[j]]<-sigma2
      }
      
      else if(j==2){
        C_sn<-sigma2*exp(-phis*neighbor.train$NN_dist[1,1])
        C_N<-sigma2
        B[[j]]<-C_sn/C_N
        F[[j]]<-sigma2- B[[j]]*C_sn
      }
      
      else{
        C_sn<-sigma2*exp(-phis*as.matrix(neighbor.train$NN_dist[j-1,1:min(j-1,M)]))
        D<-matrix(rep(0,(min(M,j-1))^2),nr=min(M,j-1),nc=min(M,j-1));
        D[lower.tri(D,diag = FALSE)]<-neighbor.train$NN_distM[j-1,1:min((j-1)*(j-2)/2,M*(M-1)/2)];
        D<-D+t(D)
        C_N<-sigma2*exp(-phis*D)#+(1e-6)*diag(dim(D)[1])
        B[[j]]<-t(C_sn)%*%solve(C_N)
        F[[j]]<-sigma2-B[[j]]%*%C_sn
      }  
    }
    
    #    build.BF<-function(j){
    
    #      if(j==1){
    #        B<-0
    #        F<-sigma2
    #      }
    
    #      else if(j==2){
    #        C_sn<-sigma2*exp(-phis*neighbor.train$NN_dist[1,1])
    #        C_N<-sigma2
    #        B<-C_sn/C_N
    #        F<-sigma2- B*C_sn
    #      }
    
    #      else{
    #        C_sn<-sigma2*exp(-phis*as.matrix(neighbor.train$NN_dist[j-1,1:min(j-1,M)]))
    #        D<-matrix(rep(0,(min(M,j-1))^2),nr=min(M,j-1),nc=min(M,j-1));
    #        D[lower.tri(D,diag = FALSE)]<-neighbor.train$NN_distM[j-1,1:min((j-1)*(j-2)/2,M*(M-1)/2)];
    #        D<-D+t(D)
    #        C_N<-sigma2*exp(-phis*D)#+(1e-6)*diag(dim(D)[1])
    #        B<-t(C_sn)%*%solve(C_N)
    #        F<-sigma2-B%*%C_sn
    #      } 
    #      return(list(B=B,F=F))
    #    }
    
    #    BF<-mclapply(1:dim(NN_ind)[1], build.BF, mc.cores=4)
    
    K<-dim(NN_ind)[1]
    construct.B<-function(i){
      if(i>1){
        Bs.star<-rep(0,K)
        Bs.star[i]<-1
        Bs.star[NN_ind[i,]]<- -B[[i]]
      }else{
        Bs.star<-c(1,rep(0,K-1))
      }
      return(Bs.star)
    }
    
    
    Bs<-mclapply(1:K,construct.B,mc.cores = 4)
    Bs<-do.call(cbind,Bs)
    #    F.inv<-bdiag(mclapply(1:K,function(i){ 1/BF[[i]]$F}))
    Q<-Matrix(Bs,sparse = T)%*%Matrix(diag(1/unlist(F)),sparse = T)%*%Matrix(t(Bs),sparse = T)
    
    return (Q)
    #    return(Q)
  }
  
  mu = function() {
    return(numeric(0))
  }
  
  log.norm.const = function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    require(invgamma)
    param = interpret.theta()
    
    res <- dinvgamma(param$sigma2, 2, 1, log = TRUE) + log(param$sigma2) +
      #dgamma(param$phis, 1, 5e-05, log = TRUE)  + log(param$phis) 
      log(b-a)+theta[2]-2*log(1+exp(theta[2]))-log(b-a)
    
    return(res)
  }
  
  initial <- function() {
    return(c(0, -3))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  if (!length(theta)) theta = initial()
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

######################################################################
library("INLA")

inla.setOption( num.threads = 4 )

a<-0.03
b<-100

nngp.model <- inla.rgeneric.define(inla.rgeneric.nngp.model, neighbor.train=neighbor.train,M=M,a=a,b=b,debug = T)

idx<-match(train$location,train.coord$location)

train$area<-as.vector(scale(train$square))


f.nngp<-log(totalPrice)~ -1 + intercept + area + livingRoom + drawingRoom + kitchen + bathRoom +
  floor + renovationCondition + elevator + subway + f(idx,model = nngp.model)

tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
m.nngp<-inla(f.nngp,data = data.frame(intercept=1,train),control.compute = list(dic=TRUE,waic=TRUE,config=TRUE))
sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

summary(m.nngp)

marg.tau<-inla.tmarginal(function(x){1/x},m.nngp$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.tau)
marg.sig<-inla.tmarginal(function(x){exp(-x)},m.nngp$marginals.hyperpar$`Theta1 for idx`)
inla.zmarginal(marg.sig)
marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.nngp$marginals.hyperpar$`Theta2 for idx`)
inla.zmarginal(marg.phis)

require(ggmap)
beijing<-get_map("beijing, China", zoom=10,maptype="hybrid")
#ggplot(data = train.coord,aes(x=Lng,y=Lat,colour=m.nngp$summary.random$idx$mean))+geom_point()+scale_color_gradient(low="black", high="red")

nngp.plot<-ggmap(beijing)+geom_point(data = train.coord,aes(x=Lng,y=Lat,colour=m.nngp$summary.random$idx$mean))+
  scale_color_gradient(name = "spatial\nresidual",low="white", high="red",limits=c(-0.84,1.27))+labs(x="longitude",y="latitude")


require(fields)
require(akima)
rh<-train.coord$Lng
rv<-train.coord$Lat
rz<-m.nngp$summary.random$idx$mean
rsp<-interp(rh,rv,rz,nx=500,ny=500,extrap=T)
image.plot(rsp)

n.sample <- 1e4
set.seed(5)
inla.pos<-inla.posterior.sample(n=n.sample,m.nngp,num.threads = 4,intern = FALSE,selection = list(intercept=1,area=1,livingRoom=1,
                                                                                                  drawingRoom=1, kitchen=1, bathRoom=1, floor=1,
                                                                                                  renovationCondition=1, elevator=1, subway=1,
                                                                                                  idx=0))

#####predict spatial effect 

pre_Z.inla<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.inla)<-test.coord$location

get.z<-function(i="ith sample",l="lth location"){ 
      distance<-spDistsN1(as.matrix(train.coord[,1:2]),as.matrix(test.coord[l,1:2]),longlat = T)
      near.distance<-sort(distance,decreasing = FALSE,index.return=TRUE)
      dist.NN<-near.distance$x[1:M] #dist to nearest neighbors
      ind.NN<-near.distance$ix[1:M] #index of nearest neighbors
      
#      distM.NN<-spDists(as.matrix(train.coord[ind.NN,1:2]),longlat = T)

      theta1<-inla.pos[[i]]$hyperpar[2]
      theta2<-inla.pos[[i]]$hyperpar[3]
      
      sigma2<-exp(-theta1)
      phis<-(a+b*exp(theta2))/(1+exp(theta2))
      
      C_sn<-sigma2*exp(-as.matrix(dist.NN)*phis)
      D<-spDists(as.matrix(train.coord[ind.NN,1:2]),longlat = T)
      C_N<-sigma2*exp(-D*phis)
      B.new<-t(C_sn)%*%solve(C_N)
      F.new<-sqrt(sigma2-B.new%*%C_sn)
      
      if(sigma2-B.new%*%C_sn<0) print("error: sigma2-B.new%*%C_sn<0")
      
      prediction.z<-rnorm(1,B.new%*%inla.pos[[i]]$latent[ind.NN],F.new)
      return(prediction.z)
    }

for (l in 1:dim(test.coord)[1]) {
  print(l)
  pos.sample<- mclapply(1:n.sample,get.z, l, mc.cores=4)
  pre_Z.inla[,l]<-unlist(pos.sample)
}

################Spatial Test############
require(akima)
require(fields)
rh<-test.coord$Lng
rv<-test.coord$Lat
rz<-colMeans(pre_Z.inla)
rsp<-interp(rh,rv,rz,nx=200,ny=200,extrap=T)
image.plot(rsp)
########################################

X <- as.matrix(cbind(1,train[,c("area","livingRoom","drawingRoom","kitchen","bathRoom","floor","renovationCondition","elevator","subway")]))
y<-log(train$totalPrice)
nzip<-dim(train.coord)[1]

get.est<-function(i){
  precision<-inla.pos[[i]]$hyperpar[1]
  res1<-rnorm(n=dim(X)[1],X%*%as.vector(inla.pos[[i]]$latent[-(1:nzip)])+rep(inla.pos[[i]]$latent[1:nzip],train.coord$Freq),sqrt(1/precision))
  return(res1)
}

y.rep<-mclapply(1:n.sample,get.est,mc.cores = 4)
y.rep<-do.call(rbind,y.rep)

test$area<-as.vector(scale(test$square,center = mean(train$square), scale = sd(train$square)))
X.test <- as.matrix(cbind(1,test[,c("area","livingRoom","drawingRoom","kitchen","bathRoom","floor","renovationCondition","elevator","subway")]))
y.test<-log(test$totalPrice)
nzip<-dim(train.coord)[1]

get.pre<-function(i){
  precision<-inla.pos[[i]]$hyperpar[1]
  res2<-rnorm(n=dim(X.test)[1],X.test%*%as.vector(inla.pos[[i]]$latent[-(1:nzip)])+rep(pre_Z.inla[i,],test.coord$Freq),sqrt(1/precision))
  return(res2)
}

y.pre<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(rbind, y.pre)

pre.inla<-colMeans(y.pre)

#inla.pos.est<-inla.posterior.sample(n=15000,m.nngp,num.threads = 4,intern = FALSE,selection = list(Predictor=0))
#y.rep<-t(sapply(inla.pos.est,function(x) x$latent))
est.inla<-colMeans(y.rep)


#RMSPE<-sqrt(mean(((pre.inla-y.test)/y.test)^2))

RMSE<- sqrt(mean((pre.inla-y.test)^2))


G<-sum((y-est.inla)^2)
P<-sum(apply(y.rep, 2, var))
D<-G+P

#######################################
##############Non-spatial##############
#######################################

f.iid<-log(totalPrice)~ -1 + intercept + area + livingRoom + drawingRoom + kitchen + bathRoom +
  floor + renovationCondition + elevator + subway + f(idx,model = "iid")

tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
m.iid<-inla(f.iid,data = data.frame(intercept=1,train),control.compute = list(dic=TRUE,waic=TRUE,config=TRUE))
sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

summary(m.iid)

set.seed(5)
n.sample <- 1e4
inla.pos.iid<-inla.posterior.sample(n=n.sample,m.iid,num.threads = 4,intern = FALSE,selection = list(intercept=1,area=1,livingRoom=1,
                                                                                                  drawingRoom=1, kitchen=1, bathRoom=1, floor=1,
                                                                                                  renovationCondition=1, elevator=1, subway=1,
                                                                                                  idx=0))


pre_Z.iid<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.iid)<-test.coord$location

get.z<-function(i="ith sample",l="lth location"){ 

  theta1<-inla.pos.iid[[i]]$hyperpar[2]

  sigma2<-1/theta1
  
  prediction.z<-rnorm(1,0,sqrt(sigma2))
  return(prediction.z)
}

for (l in 1:dim(test.coord)[1]) {
  print(l)
  pos.sample<- mclapply(1:n.sample,get.z, l, mc.cores=4)
  pre_Z.iid[,l]<-unlist(pos.sample)
}


X <- as.matrix(cbind(1,train[,c("area","livingRoom","drawingRoom","kitchen","bathRoom","floor","renovationCondition","elevator","subway")]))
y<-log(train$totalPrice)
nzip<-dim(train.coord)[1]
get.est<-function(i){
  precision<-inla.pos.iid[[i]]$hyperpar[1]
  res1<-rnorm(n=dim(X)[1],X%*%as.vector(inla.pos.iid[[i]]$latent[-(1:nzip)])+rep(inla.pos.iid[[i]]$latent[1:nzip],train.coord$Freq),sqrt(1/precision))
  return(res1)
}

y.rep<-mclapply(1:n.sample,get.est,mc.cores = 4)
y.rep<-do.call(rbind,y.rep)

test$area<-as.vector(scale(test$square,center = mean(train$square), scale = sd(train$square)))
X.test <- as.matrix(cbind(1,test[,c("area","livingRoom","drawingRoom","kitchen","bathRoom","floor","renovationCondition","elevator","subway")]))
y.test<-log(test$totalPrice)
nzip<-dim(train.coord)[1]

get.pre<-function(i){
  precision<-inla.pos.iid[[i]]$hyperpar[1]
  res2<-rnorm(n=dim(X.test)[1],X.test%*%as.vector(inla.pos.iid[[i]]$latent[-(1:nzip)])+rep(pre_Z.iid[i,],test.coord$Freq),sqrt(1/precision))
  return(res2)
}

y.pre<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(rbind, y.pre)

pre.iid<-colMeans(y.pre)

est.iid<-colMeans(y.rep)


#RMSPE<-sqrt(mean(((pre.inla-y.test)/y.test)^2))

RMSE<- sqrt(mean((pre.iid-y.test)^2))


G<-sum((y-est.iid)^2)
P<-sum(apply(y.rep, 2, var))
D<-G+P

##############################################
##############Block-NNGP######################
##############################################

##Blocking
L=7 #2^L blocks

train.coord$block.ind<-1
train.coord$ind<-1
for (i in 1:L) {
  if(i%%2==1){
    d<- 1
  }else{d<- 2}
  K<- 2^(i-1)
  for (k in 1:K) {
    cut.point<-median(train.coord[train.coord$ind==k ,d])
    train.coord[train.coord$ind==k &train.coord[,d]<=cut.point,"block.ind"]<- 2*k-1
    train.coord[train.coord$ind==k &train.coord[,d]>cut.point,"block.ind"]<-2*k
  }
  train.coord$ind<- train.coord$block.ind
}
train.coord$ind<-NULL
ggplot(data = train.coord,mapping = aes(x=Lng,y=Lat,color=as.factor(block.ind)))+geom_point()+
  geom_text(aes(label=as.factor(block.ind)),hjust=0, vjust=0)+ theme(legend.position = "none")


block.center<-cbind(aggregate(train.coord$Lng,list(train.coord$block.ind),mean)[,2],
                    aggregate(train.coord$Lat,list(train.coord$block.ind),mean)[,2],
                    1:2^L)
block.center<-as.data.frame(block.center)

M=4 #number of neighbors
NN_ind<-matrix(rep(0, (dim(block.center)[1]-1)*M),nrow = dim(block.center)[1]-1,ncol = M)
NN_distM<-matrix(rep(0, (dim(block.center)[1]-1)*M*(M-1)/2),nrow = dim(block.center)[1]-1,ncol = M*(M-1)/2)
NN_dist<-matrix(rep(0, (dim(block.center)[1]-1)*M),nrow = dim(block.center)[1]-1,ncol = M)

for (i in 2:dim(block.center)[1]) {
  distance<-spDistsN1(as.matrix(block.center[1:(i-1),1:2]),as.matrix(block.center[i,1:2]),longlat = T)
  near.distance<-sort(distance,decreasing = FALSE,index.return=TRUE)
  dist.NN<-near.distance$x[1:min(i-1,M)] #dist to nearest neighbors
  ind.NN<-near.distance$ix[1:min(i-1,M)] #index of nearest neighbors
  distM.NN<-spDists(as.matrix(block.center[ind.NN,1:2]),longlat = T)
  
  NN_ind[i-1,1:min(i-1,M)]<-ind.NN
  if(i>=3){
    NN_distM[i-1,1:(min(i-1,M)*(min(i-1,M)-1)/2)]<-distM.NN[lower.tri(distM.NN,diag = FALSE)]
  }
  NN_dist[i-1,1:min(i-1,M)]<-dist.NN
}

neighbor.block<-list(NN_ind=NN_ind, NN_distM=NN_distM, NN_dist=NN_dist)


#train.coord<-train.coord[order(match(train.coord$block.ind,neighbor.block$ord)),]
#train.coord$block.ind<-rep(1:nb^2,table(train.coord$block.ind)[neighbor.block$ord])
#ggplot(data = train.coord,mapping = aes(x=horizontal,y=vertical,color=as.factor(block.ind)))+geom_point()

train.coord<-train.coord[order(train.coord$block.ind),]
train<-train[order(match(train$location,train.coord$location)),]

D_bN<-list()
D_NN<-list()
D_b<-list()
for (i in 1:2^L) {
  if(i==1){
    D_b[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),longlat = T)
    
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
    D_bN[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),as.matrix(nn.block[, c("Lng","Lat")]),longlat = T)
    D_NN[[i]]<-spDists(as.matrix(nn.block[, c("Lng","Lat")]),longlat = T)
    D_b[[i]]<-spDists(as.matrix(train.coord[train.coord$block.ind==i, c("Lng","Lat")]),longlat = T)
  }
}


inla.rgeneric.blocknngp.model <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL) {
  
  envir = parent.env(environment())
  #Internal function
  
  interpret.theta <- function() {
    return(
      list(sigma2 = exp(-theta[1]),
           phis = (a+b*exp(theta[2]))/(1+exp(theta[2])))
    )
  }
  
  graph = function() {
    require(Matrix)
    
    return (Q())
  }
  
  Q = function() {
    require(Matrix)
    require(parallel)
    
    param = interpret.theta()
    
    sigma2<-param$sigma2
    phis<-param$phis
    
    #BLOCK LEVEL DISTRIBUTION
    build.BF<-function(i){
      if(i==1){
        B<-0
        F<-sigma2*exp(-phis*D_b[[i]])
      }else{
        C_bN<-sigma2*exp(-phis*D_bN[[i]])
        C_N<-sigma2*exp(-phis*D_NN[[i]])
        C_b<-sigma2*exp(-phis*D_b[[i]])
        B<-C_bN%*%solve(C_N)
        F<-C_b-B%*%t(C_bN)
      }
      return(list(B=B,F=F))
    }
    ####Construct Precision Matrix
    res.BF<-mclapply(1:2^nb,build.BF,mc.cores = 4)
    K<-nzip
    construct.Bs<-function(i){
      if(i>1){
        nb.i<-dim(res.BF[[i]]$F)[1]
        Bs.star<-matrix(rep(0,K*nb.i),nrow = nb.i,ncol = K)
        Bs.star[,train.coord$block.ind==i]<-diag(nb.i)
        Bs.star[,which(train.coord$block.ind %in% neighbor.block$NN_ind[i-1,])]<- -res.BF[[i]]$B
      }else{
        nb.i<-dim(res.BF[[i]]$F)[1]
        Bs.star<-matrix(rep(0,K*nb.i),nrow = nb.i,ncol = K)
        Bs.star[,train.coord$block.ind==i]<-diag(nb.i)
      }
      return(Matrix(t(Bs.star),sparse = T))
    }
    
    Bs<-mclapply(1:2^nb,construct.Bs,mc.cores = 4)
    Bs<-do.call(cbind,Bs)
    
    F.inv<-bdiag(mclapply(1:2^nb,function(i){ solve(res.BF[[i]]$F)}))
    
    Q<-Bs%*%F.inv%*%t(Bs)
    
    return (Q)
    #    return(Q)
  }
  
  mu = function() {
    return(numeric(0))
  }
  
  log.norm.const = function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    require(invgamma)
    param = interpret.theta()
    
    res <- dinvgamma(param$sigma2, 2, 1, log = TRUE) + log(param$sigma2) +
      #dgamma(param$phis, 1, 5e-05, log = TRUE)  + log(param$phis) 
      log(b-a)+theta[2]-2*log(1+exp(theta[2]))-log(b-a)
    
    return(res)
  }
  
  initial <- function() {
    return(c(0, -3))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  if (!length(theta)) theta = initial()
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

library("INLA")

inla.setOption( num.threads = 4 )

a<-0.03
b<-100

block.nngp.model <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                         nb=L, D_bN=D_bN, D_NN=D_NN, D_b=D_b, a=a, b=b,debug = T)

idx<-match(train$location,train.coord$location)

train$area<-as.vector(scale(train$square))


f.blocknngp<-log(totalPrice)~ -1 + intercept + area + livingRoom + drawingRoom + kitchen + bathRoom +
  floor + renovationCondition + elevator + subway + f(idx,model = block.nngp.model)

tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
m.blocknngp<-inla(f.blocknngp,data = data.frame(intercept=1,train),control.compute = list(dic=TRUE,waic=TRUE,config=TRUE))
sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

summary(m.blocknngp)


marg.tau<-inla.tmarginal(function(x){1/x},m.blocknngp$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.tau)
marg.sig<-inla.tmarginal(function(x){exp(-x)},m.blocknngp$marginals.hyperpar$`Theta1 for idx`)
inla.zmarginal(marg.sig)
marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.blocknngp$marginals.hyperpar$`Theta2 for idx`)
inla.zmarginal(marg.phis)

blocknngp.plot<-ggmap(beijing)+geom_point(data = train.coord,aes(x=Lng,y=Lat,colour=m.blocknngp$summary.random$idx$mean))+
  scale_color_gradient(name = "spatial\nresidual",low="white", high="red",limits=c(-0.84,1.27))+labs(x="longitude",y="latitude")




require(fields)
require(akima)
rh<-train.coord$Lng
rv<-train.coord$Lat
rz<-m.blocknngp$summary.random$idx$mean
rsp<-interp(rh,rv,rz,nx=500,ny=500,extrap=T)
image.plot(rsp)

n.sample <- 1e4
set.seed(5)
block.nngp.pos<-inla.posterior.sample(n=n.sample,m.blocknngp,num.threads = 4,intern = FALSE,selection = list(intercept=1,area=1,livingRoom=1,
                                                                                                  drawingRoom=1, kitchen=1, bathRoom=1, floor=1,
                                                                                                  renovationCondition=1, elevator=1, subway=1,
                                                                                                  idx=0))


predict.spatial<- function(j="iterations",i="ith location"){
  
  theta1<-block.nngp.pos[[j]]$hyperpar[2]
  theta2<-block.nngp.pos[[j]]$hyperpar[3]
  
  sigma2<-exp(-theta1)
  phis<-(a+b*exp(theta2))/(1+exp(theta2))
  
  b.i<-which.min(spDistsN1(as.matrix(block.center[,1:2]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
  C_iN<-sigma2*exp(-phis*spDistsN1(as.matrix(train.coord[train.coord$block.ind==b.i,c("Lng","Lat")]), as.matrix(test.coord[i,c("Lng","Lat")]),longlat = TRUE))
  C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
  B.i<-C_iN %*% C_N.inv
  F.i<-sigma2-B.i%*%C_iN
  z.i<-rnorm(1, B.i%*% block.nngp.pos[[j]]$latent[which(train.coord$block.ind==b.i)],sqrt(F.i))
  return(z.i)
}

pre_Z.block.nngp<-matrix(NA, nr=n.sample,nc=dim(test.coord)[1])
colnames(pre_Z.block.nngp)<-test.coord$location
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res1<-mclapply(1:n.sample,predict.spatial,i=i,mc.cores = 4)
  pre_Z.block.nngp[,i]<-unlist(res1)
}


################Spatial Test############
require(akima)
require(fields)
rh<-test.coord$Lng
rv<-test.coord$Lat
rz<-colMeans(pre_Z.block.nngp)
rsp<-interp(rh,rv,rz,nx=200,ny=200,extrap=T)
image.plot(rsp)
########################################



X <- as.matrix(cbind(1,train[,c("area","livingRoom","drawingRoom","kitchen","bathRoom","floor","renovationCondition","elevator","subway")]))
y<-log(train$totalPrice)
nzip<-dim(train.coord)[1]

get.est<-function(i){
  
  precision<-block.nngp.pos[[i]]$hyperpar[1]
  res1<-rnorm(n=dim(X)[1],X%*%as.vector(block.nngp.pos[[i]]$latent[-(1:nzip)])+rep(block.nngp.pos[[i]]$latent[1:nzip],train.coord$Freq),sqrt(1/precision))
  return(res1)
  
}

y.rep<-mclapply(1:n.sample,get.est,mc.cores = 4)
y.rep<-do.call(rbind,y.rep)


test$area<-as.vector(scale(test$square,center = mean(train$square), scale = sd(train$square)))
X.test <- as.matrix(cbind(1,test[,c("area","livingRoom","drawingRoom","kitchen","bathRoom","floor","renovationCondition","elevator","subway")]))
y.test<-log(test$totalPrice)
nzip<-dim(train.coord)[1]

get.pre<-function(i){
  
  precision<-block.nngp.pos[[i]]$hyperpar[1]
  res2<-rnorm(n=dim(X.test)[1],X.test%*%as.vector(block.nngp.pos[[i]]$latent[-(1:nzip)])+rep(pre_Z.block.nngp[i,],test.coord$Freq),sqrt(1/precision))
  return(res2)
}

y.pre<-mclapply(1:n.sample, get.pre, mc.cores = 4)
y.pre<-do.call(rbind, y.pre)

pre.block.nngp<-colMeans(y.pre)

#inla.pos.est<-inla.posterior.sample(n=15000,m.nngp,num.threads = 4,intern = FALSE,selection = list(Predictor=0))
#y.rep<-t(sapply(inla.pos.est,function(x) x$latent))
est.block.nngp<-colMeans(y.rep)


#RMSPE<-sqrt(mean(((pre.block.nngp-y.test)/y.test)^2))

RMSE<- sqrt(mean((pre.block.nngp-y.test)^2))


G<-sum((y-est.block.nngp)^2)
P<-sum(apply(y.rep, 2, var))
D<-G+P




