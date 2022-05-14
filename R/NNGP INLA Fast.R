
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

nngp.model <- inla.rgeneric.define(inla.rgeneric.nngp.model, neighbor.train=neighbor.train,M=10,a=a,b=b,debug = T)

idx<-match(train$zip,train.coord$zip)

f.nngp<-nclaims~ -1 + intercept + cov1 + cov2 + f(idx,model = nngp.model)

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

inla.pos<-inla.posterior.sample(n=2000,m.nngp,num.threads = 4,intern = FALSE,selection = list(intercept=1,cov1=1,cov2=1,idx=0))

#####predict spatial effect 
set.seed(5)

M=10
pre_Z.inla<-matrix(NA, nr=2000,nc=nzip.test)
colnames(pre_Z.inla)<-test.coord$zip
t1<-system.time(
  for (i in 1:2000) {
    
    for (l in 1:nzip.test) {
      
      aa<-as.matrix(pdist(test.coord[l,c('horizontal','vertical')],train.coord[,c('horizontal','vertical')]))
      bb<-sort(aa,decreasing = FALSE,index.return=TRUE)
      dist.NN<-bb$x[1:M] #dist to nearest neighbors
      ind.NN<-bb$ix[1:M] #ind of nearest neighbors
      
      theta1<-inla.pos[[i]]$hyperpar[2]
      theta2<-inla.pos[[i]]$hyperpar[3]
      
      sigma2<-exp(-theta1)
      phis<-(a+b*exp(theta2))/(1+exp(theta2))
      
      C_sn<-sigma2*exp(-as.matrix(dist.NN)*phis)
      D<-as.matrix(dist(train.coord[ind.NN,c('horizontal','vertical')]))
      C_N<-sigma2*exp(-D*phis)#+(1e-6)*diag(dim(D)[1])
      B.new<-t(C_sn)%*%solve(C_N)
      F.new<-sqrt(sigma2-B.new%*%C_sn)
      
      if(sigma2-B.new%*%C_sn<0) break
      
      prediction.z<-rnorm(1,B.new%*%inla.pos[[i]]$latent[ind.NN],F.new)
      pre_Z.inla[i,l]<-prediction.z
      
    }
  }
)

#pre<-rep(0,dim(X.test)[1])
#y.pre<-NULL
#y.rep<-NULL
#for (i in 1:5000) {
#  precision<-inla.pos[[i]]$hyperpar[1]
  
#  y.pre<-rbind(y.pre,rnorm(n=dim(X.test)[1],X.test%*%as.vector(inla.pos[[i]]$latent[(nzip+1):(nzip+3)])+rep(pre_Z.inla[i,],test.coord$nhouse),sqrt(1/precision)))
#  y.rep<-rbind(y.rep,rnorm(n=dim(X)[1],X%*%as.vector(inla.pos[[i]]$latent[(nzip+1):(nzip+3)])+rep(inla.pos[[i]]$latent[1:nzip],train.coord$nhouse),sqrt(1/precision)))
  
#}


X <- cbind(1,train$cov1,train$cov2)
y<-train$nclaims
get.est<-function(i){
  
  precision<-inla.pos[[i]]$hyperpar[1]
  res1<-rnorm(n=dim(X)[1],X%*%as.vector(inla.pos[[i]]$latent[(nzip+1):(nzip+3)])+rep(inla.pos[[i]]$latent[1:nzip],train.coord$nhouse),sqrt(1/precision))
  return(res1)
  
}

y.rep<-mclapply(1:2000,get.est,mc.cores = 4)
y.rep<-do.call(rbind,y.rep)

get.pre<-function(i){
  
  precision<-inla.pos[[i]]$hyperpar[1]
  res2<-rnorm(n=dim(X.test)[1],X.test%*%as.vector(inla.pos[[i]]$latent[(nzip+1):(nzip+3)])+rep(pre_Z.inla[i,],test.coord$nhouse),sqrt(1/precision))
  return(res2)
}

y.pre<-mclapply(1:2000, get.pre, mc.cores = 4)
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



####
#### 607*524

rh<-train.coord$horizontal
rv<-train.coord$vertical
rz<-m.nngp$summary.random$idx$mean
rsp<-interp(rh,rv,rz,nx=200,ny=200,extrap=T)
image.plot(rsp,zlim=c(-4.2,2.7))


