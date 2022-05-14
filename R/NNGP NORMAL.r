




require(Matrix)

require(parallel)

require(sp)

require(pdist)
require(invgamma)

require(akima)
require(fields)
# True parameters

beta0_true <- 2;beta1_true <- 1;beta2_true<--1

phis_true <- 10; sig2 <- 1; p.cov <- 0.5; tau<-0.1

nzip <- 2000

coord<-cbind(c(10001:(nzip+10000)),runif(nzip),runif(nzip))

nhouse<-sample(1:3,nzip,replace = T)
#nhouse<-rep(1,nzip)
coord<-cbind(coord,nhouse)
colnames(coord)<-c('zip','horizontal','vertical','nhouse')
coord<-as.data.frame(coord)
#nhouse[nhouse<=5]<-5

#nhouse[1:(N-sum(nhouse))] <- nhouse[1:(N-sum(nhouse))]+1 

distM <- as.matrix(dist(coord[,c('horizontal','vertical')])) #Euclidean 
#distM<-spDists(as.matrix(ne.zip[,4:5]), longlat = TRUE, segments = FALSE, diagonal = FALSE)/100 #Great Circle distance

R <- sig2*exp(-distM*phis_true)


cov1 <- unlist(sapply(coord$nhouse, function(x) rbinom(x,1,p.cov)))

cov2<- unlist(sapply(nhouse, function(x) rnorm(x,0,1)))

#cov3 <- unlist(sapply(nhouse, function(x) rbinom(x,1,p.cov-0.1)))

#cov4<- unlist(sapply(nhouse, function(x) rnorm(x,0,1)))


trueZ<-mvtnorm::rmvnorm(1, rep(0,nzip), sigma=R, method = "chol")
trueZ<-as.vector(trueZ)
#trueZ <- MASS::mvrnorm(1, rep(0,nzip), Sigma=R) #spatial effect
summary(trueZ)

####TRUE SPATIAL PLOT###
hh<-coord$horizontal
vv<-coord$vertical
zz<-trueZ
sp<-interp(hh,vv,zz,nx=200,ny=200)
image.plot(sp)
#######################

names(trueZ) <- coord$zip
nclaims <- sapply(beta0_true+beta1_true*cov1+beta2_true*cov2+rep(trueZ,coord$nhouse), function(x) rnorm(1,mean=x,sd=sqrt(tau)),simplify = T)
sim_data<- data.frame(zip = names(nclaims), nclaims = as.numeric(nclaims), cov1=cov1,cov2=cov2)

ind.train<-sort(sample(nzip,floor(nzip*0.8)))
ind.test<-c(1:nzip)[-ind.train]


train.coord<-coord[ind.train,]
test.coord<-coord[ind.test,]


train<-sim_data[which(sim_data$zip%in%train.coord$zip),]
test<-sim_data[which(sim_data$zip%in%test.coord$zip),]


original.train.coord<-train.coord
original.test.coord<-test.coord

M=20


neighbor.train<-NNMatrix(original.train.coord[,c('horizontal','vertical')],M,1,'cb')
train.coord<-original.train.coord[neighbor.train$ord,]
train<-train[order(match(train$zip,train.coord$zip)),]


nzip<-dim(train.coord)[1]
nzip.test<-dim(test.coord)[1]



X <- cbind(1,train$cov1,train$cov2)
y <- train$nclaims

X.test<-cbind(1,test$cov1,test$cov2)
y.test<-test$nclaims

#set.seed(5)
#initial values
beta <- beta_start <-c(rnorm(1,2,1),0,0)#coefficients(lm(y~X-1)) ##rep(0,dim(X)[2])
t<-2001
z.temp<-rep(0,nzip)

phis <- phi_start <- 10
sigma2 <- sig2_start <- 1^2
tau<-0.1

c <- 1e-3
accept_mat <- c()
accept <- rep(0,2)
tuning <- rep(0.5,2)
#R <- exp(-phis*distM^2); R.inv <- chol2inv(chol(R+c*diag(nzip))); E <- eigen(R+c*diag(nzip))$values
#  system.time(chol2inv(chol(R+c*diag(nzip))))

report <- 100; 
M=10
lbeta <- length(beta)
niter <- 2e4
res_beta <- matrix(NA, nr=niter+1,nc=lbeta)
res_Z <- matrix(NA, nr=niter+1,nc=nzip*length(t))
res_sigma2 <- rep(NA,length=niter+1)
rphis <- rep(NA,length=niter+1)
res_tau<-rep(NA,length=niter+1)

pre_Z <- matrix(NA, nr=niter+1,nc=nzip.test)
colnames(pre_Z)<-test.coord$zip
colnames(res_Z)<-train.coord$zip


NN_ind<-neighbor.train$NN_ind
NN_ind<-rbind(rep(0,dim(NN_ind)[2]),NN_ind)
zip.house<-c(0,cumsum(train.coord$nhouse))

B<-list()
F<-list()
for (j in 1:nzip) {
  
  if(j==1){
    B[[j]]<-0
    F[[j]]<-sqrt(sigma2)
  }
  
  else if(j==2){
    C_sn<-sigma2*exp(-neighbor.train$NN_dist[1,1]*phis)
    C_N<-sigma2
    B[[j]]<-C_sn/C_N
    F[[j]]<-sqrt(sigma2- B[[j]]*C_sn)
  }
  
  else{
    C_sn<-sigma2*exp(-as.matrix(neighbor.train$NN_dist[j-1,1:min(j-1,M)])*phis)
    D<-matrix(rep(0,(min(M,j-1))^2),nr=min(M,j-1),nc=min(M,j-1));
    D[lower.tri(D,diag = FALSE)]<-neighbor.train$NN_distM[j-1,1:min((j-1)*(j-2)/2,M*(M-1)/2)];
    #D[upper.tri(D,diag = FALSE)]<-t(D)[lower.tri(D,diag = FALSE)]
    D<-D+t(D)
    C_N<-sigma2*exp(-D*phis)#+(1e-6)*diag(dim(D)[1])
    B[[j]]<-t(C_sn)%*%solve(C_N)
    F[[j]]<-sqrt(sigma2-B[[j]]%*%C_sn)
  }  
  
}

#z.temp<-trueZ[neighbor.train$ord]



t3<-system.time(  
  
  for(i in 2:(niter+1)){
    
    
    if(i%%report==0){
      cat("Iteration ",i,"\n")
      cat('accept',accept,"\n")
      cat('tuning',tuning,"\n")
      accept_mat <- cbind(accept_mat,accept)
      for(avec in 1:2){
        accept[avec] <- accept[avec]/100;
        accept[avec] <- max(0.1667, min(accept[avec],0.75))
        if(accept[avec]>0.5) tuning[avec] <- tuning[avec]*accept[avec]/0.5
        else if(accept[avec]<0.25) tuning[avec] <- tuning[avec]*accept[avec]/0.25
        
      }
      accept <- rep(0,2)
      
    }
    
    
    
    
    if(i %in% c(seq(100,niter,by=200))){plot(res_beta[,1],type = 'l',main = 'beta_0');abline(h=beta0_true)}
    else if(i %in% c(seq(200,niter,by=200))){
      plot(rphis,type = 'l',main = 'phis');abline(h=phis_true)
    }else if (i %in% c(seq(150,niter,by=200))){
      plot(res_Z[,100],type = 'l',main = 'spatial');
    }
    
    
    ##UPDATE BETA
    V_beta<-solve(1/10^2*diag(length(beta))+t(X)%*%X/tau)
    mu_beta<-tau*t(X)%*%(y-rep(z.temp,train.coord$nhouse))
    
    XX_inv<-solve(t(X)%*%X)
    mu_beta<-XX_inv%*%(t(X)%*%(y-rep(z.temp,train.coord$nhouse)))
    sig_beta<-tau*XX_inv
    beta<-MASS::mvrnorm(1, mu_beta, Sigma=sig_beta)
    res_beta[i,] <- beta
    
    # update tau
    tau<-rinvgamma(1,dim(train)[1]/2+2 ,rate = sum((y-X%*%beta-rep(z.temp,train.coord$nhouse))^2)/2+0.1)
    res_tau[i]<-tau
    # update Z
    for (j in 1:nzip) {
      
      if(j==1){
        in.neighbor<-which(NN_ind==j,arr.ind = T)[,1]
        in.neighbor.v<-0
        in.neighbor.mu<-0
        if(identical(in.neighbor, integer(0))){
          V_s<-1/(train.coord$nhouse[j]/tau + 1/(F[[j]])^2)
          mu_s<-sum(y[(zip.house[j]+1):zip.house[j+1]]-X[(zip.house[j]+1):zip.house[j+1],]%*%beta)/tau 
        }else{
          for (t in in.neighbor) {
            l<-which(j==NN_ind[t,])
            in.neighbor.v<-in.neighbor.v+ (B[[t]][l])^2/(F[[t]])^2
            in.neighbor.mu<- in.neighbor.mu+B[[t]][l]/(F[[t]])^2*(z.temp[t]-B[[t]]%*%z.temp[NN_ind[t,]]+B[[t]][l]*z.temp[j])
          }
          V_s<-train.coord$nhouse[j]/tau + 1/(F[[j]])^2 + in.neighbor.v
          mu_s<-sum(y[(zip.house[j]+1):zip.house[j+1]]-X[(zip.house[j]+1):zip.house[j+1],]%*%beta)/tau + in.neighbor.mu
        }
        
      }else{
        in.neighbor<-which(NN_ind==j,arr.ind = T)[,1]
        in.neighbor.v<-0
        in.neighbor.mu<-0
        if(identical(in.neighbor, integer(0))){
          V_s<-train.coord$nhouse[j]/tau + 1/(F[[j]])^2
          mu_s<-sum(y[(zip.house[j]+1):zip.house[j+1]]-X[(zip.house[j]+1):zip.house[j+1],]%*%beta)/tau + B[[j]]%*%z.temp[NN_ind[j,]]/(F[[j]])^2
        }else{
          for (t in in.neighbor) {
            l<-which(j==NN_ind[t,])
            in.neighbor.v<-in.neighbor.v+ (B[[t]][l])^2/(F[[t]])^2
            in.neighbor.mu<- in.neighbor.mu+B[[t]][l]/(F[[t]])^2*(z.temp[t]-B[[t]]%*%z.temp[NN_ind[t,]]+B[[t]][l]*z.temp[j])
          }
          V_s<-train.coord$nhouse[j]/tau + 1/(F[[j]])^2 + in.neighbor.v
          mu_s<-sum(y[(zip.house[j]+1):zip.house[j+1]]-X[(zip.house[j]+1):zip.house[j+1],]%*%beta)/tau + B[[j]]%*%z.temp[NN_ind[j,]]/F[[j]]^2 +in.neighbor.mu
          
        }
      }
      
      z.temp[j]<-rnorm(1,mu_s/V_s,sqrt(1/V_s))
      
    }
    res_Z[i,]<-z.temp  
    
    
    
    #update sigma
    sigma2.star<-exp(log(sigma2)+rnorm(1)*tuning[1])
    
    B.star<-list()
    F.star<-list()
    for (j in 1:nzip) {
      
      if(j==1){
        B.star[[j]]<-0
        F.star[[j]]<-sqrt(sigma2.star)
      }
      
      else if(j==2){
        C_sn<-sigma2.star*exp(-neighbor.train$NN_dist[1,1]*phis)
        C_N<-sigma2.star
        B.star[[j]]<-C_sn/C_N
        F.star[[j]]<-sqrt(sigma2.star- B.star[[j]]*C_sn)
      }
      
      else{
        C_sn<-sigma2.star*exp(-as.matrix(neighbor.train$NN_dist[j-1,1:min(j-1,M)])*phis)
        D<-matrix(rep(0,(min(M,j-1))^2),nr=min(M,j-1),nc=min(M,j-1));
        D[lower.tri(D,diag = FALSE)]<-neighbor.train$NN_distM[j-1,1:min((j-1)*(j-2)/2,M*(M-1)/2)];
        #D[upper.tri(D,diag = FALSE)]<-t(D)[lower.tri(D,diag = FALSE)]
        D<-D+t(D)
        C_N<-sigma2.star*exp(-D*phis)#+(1e-6)*diag(dim(D)[1])
        B.star[[j]]<-t(C_sn)%*%solve(C_N)
        F.star[[j]]<-sqrt(sigma2.star-B.star[[j]]%*%C_sn)
      }  
      
    }
    
    
    ja<-dnorm(z.temp[1],0,F.star[[1]],log=TRUE)
    jb<-dnorm(z.temp[1],0,F[[1]],log=TRUE)
    
    for (s in 2:nzip) {
      ja<-ja+dnorm(z.temp[s],B.star[[s]]%*%z.temp[NN_ind[s,]],F.star[[s]],log = TRUE)
      jb<-jb+dnorm(z.temp[s],B[[s]]%*%z.temp[NN_ind[s,]],F[[s]],log = TRUE)
      
    }
    
    ra<-ja+dinvgamma(sigma2.star,2,rate = 1,log = TRUE)
    rb<-jb+dinvgamma(sigma2,2,rate = 1,log = TRUE)
    
    if(log(runif(1))<ra-rb){
      res_sigma2[i]<-sigma2<-sigma2.star
      accept[1]<-accept[1]+1
      B<-B.star
      F<-F.star
    }else{
      res_sigma2[i]<-sigma2
    }
    
    #update phis
    phis.star<-exp(log(phis)+rnorm(1)*tuning[2])
    
    B.star<-list()
    F.star<-list()
    for (j in 1:nzip) {
      
      if(j==1){
        B.star[[j]]<-0
        F.star[[j]]<-sqrt(sigma2)
      }
      
      else if(j==2){
        C_sn<-sigma2*exp(-neighbor.train$NN_dist[1,1]*phis.star)
        C_N<-sigma2
        B.star[[j]]<-C_sn/C_N
        F.star[[j]]<-sqrt(sigma2- B.star[[j]]*C_sn)
      }
      
      else{
        C_sn<-sigma2*exp(-as.matrix(neighbor.train$NN_dist[j-1,1:min(j-1,M)])*phis.star)
        D<-matrix(rep(0,(min(M,j-1))^2),nr=min(M,j-1),nc=min(M,j-1));
        D[lower.tri(D,diag = FALSE)]<-neighbor.train$NN_distM[j-1,1:min((j-1)*(j-2)/2,M*(M-1)/2)];
        #D[upper.tri(D,diag = FALSE)]<-t(D)[lower.tri(D,diag = FALSE)]
        D<-D+t(D)
        C_N<-sigma2*exp(-D*phis.star)#+(1e-6)*diag(dim(D)[1])
        B.star[[j]]<-t(C_sn)%*%solve(C_N)
        F.star[[j]]<-sqrt(sigma2-B.star[[j]]%*%C_sn)
      }  
      
    }
    
    ja<-dnorm(z.temp[1],0,F.star[[1]],log=TRUE)
    jb<-dnorm(z.temp[1],0,F[[1]],log=TRUE)
    
    for (s in 2:nzip) {
      ja<-ja+dnorm(z.temp[s],B.star[[s]]%*%z.temp[NN_ind[s,]],F.star[[s]],log = TRUE)
      jb<-jb+dnorm(z.temp[s],B[[s]]%*%z.temp[NN_ind[s,]],F[[s]],log = TRUE)
      
    }
    
    ra<-ja+log(phis.star<30)
    rb<-jb
    
    
    if(log(runif(1))<ra-rb){
      rphis[i]<-phis<-phis.star
      accept[2]<-accept[2]+1
      B<-B.star
      F<-F.star
    }else{
      rphis[i]<-phis
    }
    
    
    if(i>=5e3){
      
      for (l in 1:nzip.test) {
        
        aa<-as.matrix(pdist(test.coord[l,c('horizontal','vertical')],train.coord[,c('horizontal','vertical')]))
        bb<-sort(aa,decreasing = FALSE,index.return=TRUE)
        dist.NN<-bb$x[1:M] #dist to nearest neighbors
        ind.NN<-bb$ix[1:M] #ind of nearest neighbors
        
        C_sn<-sigma2*exp(-as.matrix(dist.NN)*phis)
        D<-as.matrix(dist(train.coord[ind.NN,c('horizontal','vertical')]))
        C_N<-sigma2*exp(-D*phis)#+(1e-6)*diag(dim(D)[1])
        B.new<-t(C_sn)%*%solve(C_N)
        F.new<-sqrt(sigma2-B.new%*%C_sn)
        
        if(sigma2-B.new%*%C_sn<0) break
        
        prediction.z<-rnorm(1,B.new%*%z.temp[ind.NN],F.new)
        pre_Z[i,l]<-prediction.z
        
      }
    }
    
    
  }
  
)  




cor(colMeans(pre_Z[5e3:niter,]),trueZ[names(trueZ)%in%colnames(pre_Z)])


mean((colMeans(pre_Z[5e3:niter,])-trueZ[names(trueZ)%in%colnames(pre_Z)])^2)



n=1
est <-rep(0,dim(X)[1])
pre<-rep(0,dim(X.test)[1])
y.rep<-NULL
y.pre<-NULL
for (i in seq(5e3,niter)) {
  print(i)
  y.rep<-rbind(y.rep,rnorm(n=dim(X)[1],X%*%res_beta[i,]+rep(res_Z[i,],train.coord$nhouse),sqrt(res_tau[i])))
  y.pre<-rbind(y.pre,rnorm(n=dim(X.test)[1],X.test%*%res_beta[i,]+rep(pre_Z[i,],test.coord$nhouse),sqrt(res_tau[i])))
  n<-n+1
}
pre<-colMeans(y.pre)
est<-colMeans(y.rep)
#Deviance
Dev<-function(y,mu,tau){
  dev<- -2*sum(dnorm(y,mu,sqrt(tau),log = T))
  return(dev)
}

n=1
D_bar<-0
for (i in 5e3:niter) {
  D_bar<-D_bar+(Dev(y,X%*%res_beta[i,]+rep(res_Z[i,],train.coord$nhouse),res_tau[i])-D_bar)/n
  n<-n+1
}

Pd<-D_bar-Dev(y,X%*%colMeans(res_beta[5e3:niter,])+rep(colMeans(res_Z[5e3:niter,]),train.coord$nhouse),mean(res_tau[5e3:niter]))

DIC<-Pd+D_bar

RMSPE<-sqrt(mean(((pre-y.test)/y.test)^2))

RMSE<- sqrt(mean((pre-y.test)^2))


G<-sum((y-est)^2)
P<-sum(apply(y.rep, 2, var))
D<-G+P

#G<-sum((y.test-pre)^2)
#P<-sum(apply(y.rep, 2, var))
#D<-G+P
#######################################################


apply(res_beta[5e3:niter,],2,quantile,c(0.025,0.5,0.975))
quantile(rphis[5e3:niter],c(0.025,0.5,0.975))
quantile(res_sigma2[5e3:niter],c(0.025,0.5,0.975))
quantile(res_tau[5e3:niter],c(0.025,0.5,0.975))

plot(res_beta[5e3:niter,1],type="l")
plot(res_beta[5e3:niter,2],type="l")
plot(res_beta[5e3:niter,3],type="l")

plot(rphis[5e3:niter],type="l")
plot(res_sigma2[5e3:niter],type="l")
plot(res_tau[5e3:niter],type="l")


####TEST SPATIAL PLOT###
rh<-test.coord$horizontal
rv<-test.coord$vertical
rz<-trueZ[names(trueZ)%in%colnames(pre_Z)]
rsp<-interp(rh,rv,rz,nx=200,ny=200)
image.plot(rsp, zlim=range(-3.2,3))

th<-test.coord$horizontal
tv<-test.coord$vertical
tz<-colMeans(pre_Z[5e3:niter,])
tsp<-interp(th,tv,tz,nx=200,ny=200)
image.plot(tsp, zlim=range(-3.2,3))


#####TRAIN SPATIAL PLOT#### 607*524

train_Z<-trueZ[names(trueZ) %in% train.coord$zip]
train_Z<-train_Z[order(match(names(train_Z),train.coord$zip))]



rh<-train.coord$horizontal
rv<-train.coord$vertical
rz<-train_Z
rsp<-interp(rh,rv,rz,nx=200,ny=200)
image.plot(rsp,zlim=range(-4.2,2.7))


th<-train.coord$horizontal
tv<-train.coord$vertical
tz<-colMeans(res_Z[5e3:niter,])
tsp<-interp(th,tv,tz,nx=200,ny=200)
image.plot(tsp,zlim=range(-3.1,3.1))

par(mfrow=c(1,2)) 
image.plot(rsp,zlim=range(-3,3))
image.plot(tsp,zlim=range(-3,3))
par(mfrow=c(1,1)) 


#save(coord,trueZ,sim_data,nzip,nzip.test,ind.train,ind.test,
#     file = "bml sim.RData")




