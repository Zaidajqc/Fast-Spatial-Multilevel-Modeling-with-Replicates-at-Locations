###########SPLIT BLOCKS#####
require(ggplot2)
nb<-12
v_cut <- as.numeric(cut(train.coord$vertical, nb))
h_cut <- as.numeric(cut(train.coord$horizontal, nb))

block.ind<-interaction(v_cut,h_cut)
levels(block.ind)<-1:nb^2
train.coord <- cbind(train.coord,block.ind)
train.coord<-train.coord[order(train.coord$block.ind),]
dist_0<-aggregate(train.coord$horizontal^2+train.coord$vertical^2,list(train.coord$block.ind),mean)
train.coord$block.ind<-as.factor(rep(match(dist_0$x,sort(dist_0$x)),table(train.coord$block.ind)))
#ggplot(data = train.coord,mapping = aes(x=horizontal,y=vertical,color=as.factor(block.ind)))+geom_point()

block.center<-cbind(1:nb^2,aggregate(train.coord$horizontal,list(train.coord$block.ind),mean)[,2],
                    aggregate(train.coord$vertical,list(train.coord$block.ind),mean)[,2])

M=4 ####Number of nearest blocks

neighbor.block<-NNMatrix(block.center[,2:3],M,1,'cb',ord = order(block.center[,2]+block.center[,3]))

train.coord<-train.coord[order(match(train.coord$block.ind,neighbor.block$ord)),]
train.coord$block.ind<-rep(1:nb^2,table(train.coord$block.ind)[neighbor.block$ord])
#ggplot(data = train.coord,mapping = aes(x=horizontal,y=vertical,color=as.factor(block.ind)))+geom_point()

train<-train[order(match(train$zip,train.coord$zip)),]

D_bN<-list()
D_NN<-list()
D_b<-list()
for (i in 1:nb^2) {
  if(i==1){
    D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, c("vertical","horizontal")]))
    
  }else{
    nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
    D_bN[[i]]<-as.matrix(pdist(train.coord[train.coord$block.ind==i, c("vertical","horizontal")],
                               nn.block[, c("vertical","horizontal")]))
    D_NN[[i]]<-as.matrix(dist(nn.block[,c("vertical","horizontal")]))
    D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, c("vertical","horizontal")]))
  }
}

###########BLOCK NNGP########

###FOR LOOP CONSTRUCT B & F ####SLOWER
#B<-list()
#F<-list()
#for (i in 1:nb^2) {
#  if(i==1){
#    B[[i]]<-0
#    F[[i]]<-sigma2*exp(-D_b[[i]]*phis)
#  }else{
#    C_bN<-sigma2*exp(-D_bN[[i]]*phis)
#    C_N<-sigma2*exp(-D_NN[[i]]*phis)
#    C_b<-sigma2*exp(-D_b[[i]]*phis)
#    B[[i]]<-C_bN%*%solve(C_N)
#    F[[i]]<-C_b-B[[i]]%*%t(C_bN)
#  }
#}
##############################


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
    res.BF<-mclapply(1:nb^2,build.BF,mc.cores = 4)
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
    
    Bs<-mclapply(1:nb^2,construct.Bs,mc.cores = 4)
    Bs<-do.call(cbind,Bs)
    
    F.inv<-bdiag(mclapply(1:nb^2,function(i){ solve(res.BF[[i]]$F)}))
    
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
    return(c(0, -1))
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

a<-1
b<-30

block.nngp.model <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord, M=M, neighbor.block=neighbor.block,
                                   nb=nb, D_bN=D_bN, D_NN=D_NN, D_b=D_b, a=a, b=b,debug = T)

idx<-match(train$zip,train.coord$zip)

f.nngp<-nclaims~ -1 + intercept + cov1 + cov2 + f(idx,model = block.nngp.model)

tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
m.blocknngp<-inla(f.nngp,data = data.frame(intercept=1,train),control.compute = list(dic=TRUE,waic=TRUE,config=TRUE))
sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

summary(m.blocknngp)

marg.tau<-inla.tmarginal(function(x){1/x},m.blocknngp$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(marg.tau)
marg.sig<-inla.tmarginal(function(x){exp(-x)},m.blocknngp$marginals.hyperpar$`Theta1 for idx`)
inla.zmarginal(marg.sig)
marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.blocknngp$marginals.hyperpar$`Theta2 for idx`)
inla.zmarginal(marg.phis)


rh<-train.coord$horizontal
rv<-train.coord$vertical
rz<-m.blocknngp$summary.random$idx$mean
rsp<-interp(rh,rv,rz,nx=200,ny=200)
image.plot(rsp,zlim=range(-4.2,2.7))




block.nngp.pos<-inla.posterior.sample(n=2000,m.blocknngp,num.threads = 4,intern = FALSE,selection = list(intercept=1,cov1=1,cov2=1,idx=0))


block.center<-block.center[neighbor.block$ord,]
block.center[,1]<-1:nb^2
#####predict spatial effect 
set.seed(5)

#######Locations#####
ggplot(data = train.coord,mapping = aes(x=horizontal,y=vertical,color=as.factor(block.ind)))+geom_point()+
  geom_text(aes(label=as.factor(block.ind)),hjust=0, vjust=0)+ theme(legend.position = "none")

ggplot(mapping = aes(block.center[,2],block.center[,3],color=as.factor(block.center[,1])))+geom_point()+
  geom_text(aes(label=as.factor(block.center[,1])),hjust=0, vjust=0)+ theme(legend.position = "none")
######################



predict.spatial<- function(j="iterations",i="ith location"){
  
theta1<-block.nngp.pos[[j]]$hyperpar[2]
theta2<-block.nngp.pos[[j]]$hyperpar[3]

sigma2<-exp(-theta1)
phis<-(a+b*exp(theta2))/(1+exp(theta2))

b.i<-which.min(as.matrix(pdist(test.coord[i,c("horizontal","vertical")],block.center[,2:3])))
C_iN<-sigma2*exp(-phis*as.matrix(pdist(test.coord[i,c("horizontal","vertical")],train.coord[train.coord$block.ind==b.i,c("horizontal","vertical")])))
C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
B.i<-C_iN %*% C_N.inv
F.i<-sigma2-B.i%*%t(C_iN)
z.i<-rnorm(1, B.i%*% block.nngp.pos[[j]]$latent[which(train.coord$block.ind==b.i)],sqrt(F.i))
return(z.i)
}

pre_Z.block.nngp<-matrix(NA, nr=2000,nc=nzip.test)
colnames(pre_Z.block.nngp)<-test.coord$zip
for (i in 1:dim(test.coord)[1]) {
  print(i)
  res1<-mclapply(1:2000,predict.spatial,i=i,mc.cores = 3)
  pre_Z.block.nngp[,i]<-unlist(res1)
}

X <- cbind(1,train$cov1,train$cov2)
y<-train$nclaims
get.est<-function(i){
  
  precision<-block.nngp.pos[[i]]$hyperpar[1]
  res1<-rnorm(n=dim(X)[1],X%*%as.vector(block.nngp.pos[[i]]$latent[(nzip+1):(nzip+3)])+rep(block.nngp.pos[[i]]$latent[1:nzip],train.coord$nhouse),sqrt(1/precision))
  return(res1)
  
}

y.rep<-mclapply(1:2000,get.est,mc.cores = 4)
y.rep<-do.call(rbind,y.rep)

get.pre<-function(i){
  
  precision<-block.nngp.pos[[i]]$hyperpar[1]
  res2<-rnorm(n=dim(X.test)[1],X.test%*%as.vector(block.nngp.pos[[i]]$latent[(nzip+1):(nzip+3)])+rep(pre_Z.block.nngp[i,],test.coord$nhouse),sqrt(1/precision))
  return(res2)
}

y.pre<-mclapply(1:2000, get.pre, mc.cores = 4)
y.pre<-do.call(rbind, y.pre)

pre.block.nngp<-colMeans(y.pre)

#inla.pos.est<-inla.posterior.sample(n=15000,m.nngp,num.threads = 4,intern = FALSE,selection = list(Predictor=0))
#y.rep<-t(sapply(inla.pos.est,function(x) x$latent))
est.block.nngp<-colMeans(y.rep)


RMSPE<-sqrt(mean(((pre.block.nngp-y.test)/y.test)^2))

RMSE<- sqrt(mean((pre.block.nngp-y.test)^2))


G<-sum((y-est.block.nngp)^2)
P<-sum(apply(y.rep, 2, var))
D<-G+P


