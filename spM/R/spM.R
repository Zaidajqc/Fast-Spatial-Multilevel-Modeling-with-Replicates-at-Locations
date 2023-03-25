#' Preprocess data for NNGP model
#'
#' Preprocess the data to make it meet the requirement for later inference.
#' @param data The original data.
#' @param coord n*2 coordinates corresponding to the original data.
#' @param neighbor The number of neighbor needed for NNGP model, default is 10.
#' @param order.by The ordering coordinate, value=1 means order the locations by coord[,1].
#' @return Function that returns a list.
#' \itemize{
#' \item data - Reordered data with original order index and location index.
#' \item unique.coordinates - Coordinates of unique locations.
#' \item nearest.neighbor - Information about nearest neighbors.
#' \item M - The number of user-defined neighbors.
#' }
#' @export
preprocess.nngp<-function(data, coord, neighbor=10, order.by=1){

  # coord represents the n*2 coordinates corresponding to data
  # neighbor represents the number of neighbor needed, default is 10
  # order.by represents the ordering coordinate, value=1 means order the locations by coord[,1]


#  data<-cbind(data,coord)
  data<-data.frame(index=1:dim(data)[1],data)
  data<-data[order(coord[,1],coord[,2]),]
  coord<-coord[order(coord[,1],coord[,2]),]

  require(dplyr)
  unique.coord<-rename(count(coord, coord[,1], coord[,2]), Freq = n)
  unique.coord<-data.frame(unique.coord, location=10000+(1:dim(unique.coord)[1]))
  data$location<-rep(unique.coord$location,unique.coord$Freq)

  library(spNNGP)       # Build neighbor index
  #### distance matrix for location i and its neighbors ####
  i_dist <- function(i, neighbor_index, s){
    dist(s[c(i, neighbor_index[[i - 1]]), ])
  }

  ###Neighbor distcance###
  get_NN_distM <- function (ind, ind_distM_d, M) {
    if (ind < M ){l = ind } else {l = M};
    M_i <- rep(0, M * (M - 1) / 2);
    if (l == 1) {}
    else{
      M_i[1: (l * (l - 1) / 2)] <-
        c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
    }
    return(M_i)
  }
  ###distance between current location and neighbors###
  get_NN_dist <- function (ind, ind_distM_d, M) {
    if (ind < M ){l = ind } else {l = M};
    D_i <- rep(0, M);
    D_i[1:l] <- c(ind_distM_d[[ind]])[1:l]
    return(D_i)
  }
  get_NN_ind <- function (ind, ind_distM_i, M) {
    if (ind < M ){l = ind } else {l = M};
    D_i <- rep(0, M);
    D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
    return(D_i)
  }
  NNMatrix <- function(coords, n.neighbors, n.omp.threads = 1,
                       search.type = "cb", ord = order(coords[, 1])){

    N <- nrow(coords)
    m.c <- spConjNNGP(rep(0, N) ~ 1, coords = coords,
                      n.neighbors = n.neighbors,
                      theta.alpha = c("phi" = 1, "alpha" = 1),
                      sigma.sq.IG = c(2, 1),
                      cov.model = "exponential",
                      n.omp.threads = n.omp.threads,
                      search.type = search.type,
                      ord = ord,
                      return.neighbor.info = T, fit.rep = FALSE,
                      verbose = FALSE)
    M = n.neighbors
    NN_ind <- t(sapply(1: (N - 1), get_NN_ind, m.c$neighbor.info$n.indx[-1], M))

    neighbor_dist <- sapply(2:N, i_dist, m.c$neighbor.info$n.indx[-1],
                            m.c$coords[m.c$neighbor.info$ord, ])

    NN_distM <- t(sapply(1: (N - 1), get_NN_distM, neighbor_dist, M))
    NN_dist <- t(sapply(1: (N - 1), get_NN_dist, neighbor_dist, M))

    return(list(ord = m.c$neighbor.info$ord,
                coords.ord = m.c$coords[m.c$neighbor.info$ord, ],
                NN_ind = NN_ind, NN_distM = NN_distM, NN_dist = NN_dist))
  }

  M=neighbor
  if(order.by==2){
    nearest.neighbor<-NNMatrix(coords = unique.coord[,1:2],M,1,'cb',ord= order(unique.coord[, 2]))
  }else{
    nearest.neighbor<-NNMatrix(coords = unique.coord[,1:2],M,1,'cb',ord= order(unique.coord[, 1]))
  }
  unique.coord<-unique.coord[nearest.neighbor$ord,]
  data<-data[order(match(data$location,unique.coord$location)),]

  return(list(data=data, unique.coord=unique.coord, nearest.neighbor=nearest.neighbor, M=M))

}

#' Preprocess data for block-NNGP model
#'
#' Preprocess the data to make it meet the requirement for later inference.
#' @param data The original data.
#' @param coord n*2 coordinates corresponding to the original data.
#' @param neighbor The number of neighboring blocks needed for block-NNGP model, default is 4.
#' @param cut.by 2*1 vector tell how to cut the block, e.g c(5, 7) cuts the spatial domain to 35 rectangles.
#' @return Function that returns a list.
#' \itemize{
#' \item data - Reordered data with original order index and location index.
#' \item unique.coordinates - Coordinates of unique locations.
#' \item nearest.neighbor - Information about nearest neighbors.
#' \item M - The number of user-defined neighboring blocks.
#' \item block.center - The center of partitioned blocks.
#' }
#' @export
preprocess.blocknngp<-function(data, coord, cut.block=c(6,6) ,neighbor=4){

  # coord represents the n*2 coordinates corresponding to data
  # cut.block represents how to cut the block, c(cut coord[,1], cut coord[,2])

  # neighbor represents the number of neighbor needed, default is 10
  # order.by represents the ordering coordinate, value=1 means order the locations by coord[,1]

  library(spNNGP)       # Build neighbor index
  #### distance matrix for location i and its neighbors ####
  i_dist <- function(i, neighbor_index, s){
    dist(s[c(i, neighbor_index[[i - 1]]), ])
  }

  ###Neighbor distcance###
  get_NN_distM <- function (ind, ind_distM_d, M) {
    if (ind < M ){l = ind } else {l = M};
    M_i <- rep(0, M * (M - 1) / 2);
    if (l == 1) {}
    else{
      M_i[1: (l * (l - 1) / 2)] <-
        c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
    }
    return(M_i)
  }
  ###distance between current location and neighbors###
  get_NN_dist <- function (ind, ind_distM_d, M) {
    if (ind < M ){l = ind } else {l = M};
    D_i <- rep(0, M);
    D_i[1:l] <- c(ind_distM_d[[ind]])[1:l]
    return(D_i)
  }
  get_NN_ind <- function (ind, ind_distM_i, M) {
    if (ind < M ){l = ind } else {l = M};
    D_i <- rep(0, M);
    D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
    return(D_i)
  }
  NNMatrix <- function(coords, n.neighbors, n.omp.threads = 1,
                       search.type = "cb", ord = order(coords[, 1])){

    N <- nrow(coords)
    m.c <- spConjNNGP(rep(0, N) ~ 1, coords = coords,
                      n.neighbors = n.neighbors,
                      theta.alpha = c("phi" = 1, "alpha" = 1),
                      sigma.sq.IG = c(2, 1),
                      cov.model = "exponential",
                      n.omp.threads = n.omp.threads,
                      search.type = search.type,
                      ord = ord,
                      return.neighbor.info = T, fit.rep = FALSE,
                      verbose = FALSE)
    M = n.neighbors
    NN_ind <- t(sapply(1: (N - 1), get_NN_ind, m.c$neighbor.info$n.indx[-1], M))

    neighbor_dist <- sapply(2:N, i_dist, m.c$neighbor.info$n.indx[-1],
                            m.c$coords[m.c$neighbor.info$ord, ])

    NN_distM <- t(sapply(1: (N - 1), get_NN_distM, neighbor_dist, M))
    NN_dist <- t(sapply(1: (N - 1), get_NN_dist, neighbor_dist, M))

    return(list(ord = m.c$neighbor.info$ord,
                coords.ord = m.c$coords[m.c$neighbor.info$ord, ],
                NN_ind = NN_ind, NN_distM = NN_distM, NN_dist = NN_dist))
  }


  #  data<-cbind(data,coord)
  data<-data.frame(index=1:dim(data)[1],data)
  data<-data[order(coord[,1],coord[,2]),]
  coord<-coord[order(coord[,1],coord[,2]),]

  require(dplyr)
  unique.coord<-rename(count(coord, coord[,1], coord[,2]), Freq = n)
  unique.coord<-data.frame(unique.coord, location=10000+(1:dim(unique.coord)[1]))
  data$location<-rep(unique.coord$location,unique.coord$Freq)

  nb<-cut.block
  v_cut <- as.numeric(cut(unique.coord[,1], nb[1]))
  h_cut <- as.numeric(cut(unique.coord[,2], nb[2]))

  block.ind<-interaction(v_cut,h_cut)
  levels(block.ind)<-1:(nb[1]*nb[2])
  unique.coord <- cbind(unique.coord,block.ind)
  unique.coord<-unique.coord[order(unique.coord$block.ind),]
  dist_0<-aggregate(unique.coord[,1]^2+unique.coord[,2]^2,list(unique.coord$block.ind),mean)
  unique.coord$block.ind<-as.factor(rep(match(dist_0$x,sort(dist_0$x)),table(unique.coord$block.ind)))
  #ggplot(data = train.coord,mapping = aes(x=horizontal,y=vertical,color=as.factor(block.ind)))+geom_point()

  block.center<-cbind(1:(nb[1]*nb[2]),aggregate(unique.coord[,1],list(unique.coord$block.ind),mean)[,2],
                      aggregate(unique.coord[,2],list(unique.coord$block.ind),mean)[,2])

  M=neighbor ####Number of nearest blocks

  neighbor.block<-NNMatrix(block.center[,2:3],M,1,'cb',ord = order(block.center[,2]+block.center[,3]))

  unique.coord<-unique.coord[order(match(unique.coord$block.ind,neighbor.block$ord)),]
  unique.coord$block.ind<-rep(1:(nb[1]*nb[2]),table(unique.coord$block.ind)[neighbor.block$ord])
  #ggplot(data = train.coord,mapping = aes(x=horizontal,y=vertical,color=as.factor(block.ind)))+geom_point()

  data<-data[order(match(data$location,unique.coord$location)),]

  return(list(data=data, unique.coord=unique.coord, neighbor.block=neighbor.block,M=M,cut.block=nb,block.center=block.center))
}



#' Spatial multilevel model with NNGP prior
#'
#' Fit the spatial multilevel model with NNGP prior to the preprocessed data from \code{\link{preprocess.nngp}}.
#' @param formula Model formula, e.g. y~x.
#' @param preprocess.data Preprocessed data from \code{\link{preprocess.nngp}}.
#' @param family A string indicating the likelihood family. The default is gaussian with identity link.
#' See names(inla.models()$likelihood) for a list of possible alternatives.
#' @param E Known component in the mean for the Poisson likelihoods defined as Eexp(lambda)
#' @param offset This argument is used to specify an a-priori known and fixed component to be included in the linear predictor during fitting.
#' @param phi.prior 2*1 vector c(a,b) indicates the uniform prior U(a, b) for spatial decay parameter phi.
#' @param sigma2.prior 2*1 vector c(a,b) indicates the Inverse-Gamma prior IG(a, b) for spatial marginal variance parameter sigma^2.
#' @param mc.cores The number of CPU cores should be used. For Windows user, only 1 core is allowed.
#' @param control... See \code{\link[INLA]{inla}}.
#' @return Function that returns a list
#' \itemize{
#' \item res - A fitted INLA model, see \code{\link[INLA]{inla}}.
#' \item Summary statistics of hyperparameters.

#' }
#' @export
nngp.inla<-function(formula, preprocess.data, family="gaussian", E=rep(1,dim(preprocess.data$data)[1]), offset = NULL,
                    phi.prior=c(1,30),sigma2.prior=c(2,1), mc.cores=1,
                    control.compute = list(dic=TRUE,waic=TRUE,config=TRUE),
                    control.predictor = list(),
                    control.family = list(),
                    control.inla = list(),
                    control.results = list(),
                    control.fixed = list()){

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

      res <- dinvgamma(param$sigma2, sigma2.prior[1], sigma2.prior[2], log = TRUE) + log(param$sigma2) +
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


  require("INLA")

  inla.setOption( num.threads = mc.cores )

  a<-phi.prior[1]
  b<-phi.prior[2]
  sigma2.prior=sigma2.prior

  nngp.model <- inla.rgeneric.define(inla.rgeneric.nngp.model, neighbor.train=preprocess.data$nearest.neighbor ,M=preprocess.data$M,
                                     a=a,b=b,sigma2.prior=sigma2.prior,debug = F)

  idx<-match(preprocess.data$data$location, preprocess.data$unique.coord$location)

#  f.nngp<-formula + f(idx,model = nngp.model)
  f.nngp<- update(formula, ~. + f(idx,model = nngp.model))

  if(family=="poisson"){
    m.nngp<-inla(f.nngp,data = preprocess.data$data, family = family, E=E, offset = offset,
               control.compute = control.compute,
               control.predictor = control.predictor,
               control.family = control.family,
               control.inla = ontrol.inla,
               control.results = control.results,
               control.fixed = control.fixed)
  }else{
    m.nngp<-inla(f.nngp,data = preprocess.data$data, family = family, offset = offset,
                 control.compute = control.compute,
                 control.predictor = control.predictor,
                 control.family = control.family,
                 control.inla = ontrol.inla,
                 control.results = control.results,
                 control.fixed = control.fixed)
  }

  if(family=="gaussian"){
  marg.tau<-inla.tmarginal(function(x){1/x},m.nngp$marginals.hyperpar$`Precision for the Gaussian observations`)
  tau.summary<-as.matrix(inla.zmarginal(marg.tau,silent = T))
  marg.sig<-inla.tmarginal(function(x){exp(-x)},m.nngp$marginals.hyperpar$`Theta1 for idx`)
  sigma.summary<-as.matrix(inla.zmarginal(marg.sig,silent = T))
  marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.nngp$marginals.hyperpar$`Theta2 for idx`)
  phi.summary<-as.matrix(inla.zmarginal(marg.phis,silent = T))

  return(list(res=m.nngp,tau.summary=tau.summary,sigma.summary=sigma.summary,phi.summary=phi.summary))
  }else{

    marg.sig<-inla.tmarginal(function(x){exp(-x)},m.nngp$marginals.hyperpar$`Theta1 for idx`)
    sigma.summary<-as.matrix(inla.zmarginal(marg.sig,silent = T))
    marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.nngp$marginals.hyperpar$`Theta2 for idx`)
    phi.summary<-as.matrix(inla.zmarginal(marg.phis,silent = T))

    return(list(res=m.nngp,sigma.summary=sigma.summary,phi.summary=phi.summary))

  }
}






#' Spatial multilevel model with block-NNGP prior
#'
#' Fit the spatial multilevel model with block-NNGP prior to the preprocessed data from \code{\link{preprocess.blocknngp}}.
#' @param formula Model formula, e.g. y~x.
#' @param preprocess.data Preprocessed data from \code{\link{preprocess.blocknngp}}.
#' @param family A string indicating the likelihood family. The default is gaussian with identity link.
#' See names(inla.models()$likelihood) for a list of possible alternatives.
#' @param E Known component in the mean for the Poisson likelihoods defined as Eexp(lambda)
#' @param offset This argument is used to specify an a-priori known and fixed component to be included in the linear predictor during fitting.
#' @param phi.prior 2*1 vector c(a,b) indicates the uniform prior U(a, b) for spatial decay parameter phi.
#' @param sigma2.prior 2*1 vector c(a,b) indicates the Inverse-Gamma prior IG(a, b) for spatial marginal variance parameter sigma^2.
#' @param mc.cores The number of CPU cores should be used. For Windows user, only 1 core is allowed.
#' @param control... See \code{\link[INLA]{inla}}.
#' @return Function that returns a list
#' \itemize{
#' \item res - A fitted INLA model, see \code{\link[INLA]{inla}}.
#' \item Summary statistics of hyperparameters.
#' }
#' @export
blocknngp.inla<-function(formula, preprocess.data, family="gaussian", E=rep(1,dim(preprocess.data$data)[1]), offset = NULL,
                    phi.prior=c(1,30),sigma2.prior=c(2,1), mc.cores=1,
                    control.compute = list(dic=TRUE,waic=TRUE,config=TRUE),
                    control.predictor = list(),
                    control.family = list(),
                    control.inla = list(),
                    control.results = list(),
                    control.fixed = list()){

  require(pdist)
  nb<-preprocess.data$cut.block
  train.coord<-preprocess.data$unique.coord
  neighbor.block<-preprocess.data$neighbor.block
  D_bN<-list()
  D_NN<-list()
  D_b<-list()
  for (i in 1:(nb[1]*nb[2])) {
    if(i==1){
      D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, 1:2]))

    }else{
      nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
      D_bN[[i]]<-as.matrix(pdist(train.coord[train.coord$block.ind==i, 1:2],
                                 nn.block[, 1:2]))
      D_NN[[i]]<-as.matrix(dist(nn.block[,1:2]))
      D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, 1:2]))
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
      res.BF<-mclapply(1:(nb[1]*nb[2]),build.BF,mc.cores = mc.cores)
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

      Bs<-mclapply(1:(nb[1]*nb[2]),construct.Bs,mc.cores = mc.cores)
      Bs<-do.call(cbind,Bs)

      F.inv<-bdiag(mclapply(1:(nb[1]*nb[2]),function(i){ solve(res.BF[[i]]$F)}))

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

      res <- dinvgamma(param$sigma2, sigma2.prior[1], sigma2.prior[2], log = TRUE) + log(param$sigma2) +
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

  inla.setOption( num.threads = mc.cores )


  a<-phi.prior[1]
  b<-phi.prior[2]
  sigma2.prior=sigma2.prior
  M=preprocess.data$M

  block.nngp.model <- inla.rgeneric.define(inla.rgeneric.blocknngp.model, nzip=dim(train.coord)[1], train.coord=train.coord,
                                           M=M, neighbor.block=neighbor.block,
                                           nb=nb, D_bN=D_bN, D_NN=D_NN, D_b=D_b, a=a, b=b,sigma2.prior=sigma2.prior,
                                           debug = F,mc.cores=mc.cores)

  idx<-match(preprocess.data$data$location, preprocess.data$unique.coord$location)

  #  f.nngp<-formula + f(idx,model = nngp.model)
  f.blocknngp<- update(formula, ~. + f(idx,model = block.nngp.model))

  if(family=="poisson"){
    m.blocknngp<-inla(f.blocknngp,data = preprocess.data$data, family = family, E=E, offset = offset,
                 control.compute = control.compute,
                 control.predictor = control.predictor,
                 control.family = control.family,
                 control.inla = ontrol.inla,
                 control.results = control.results,
                 control.fixed = control.fixed)
  }else{
    m.blocknngp<-inla(f.blocknngp,data = preprocess.data$data, family = family, offset = offset,
                 control.compute = control.compute,
                 control.predictor = control.predictor,
                 control.family = control.family,
                 control.inla = ontrol.inla,
                 control.results = control.results,
                 control.fixed = control.fixed)
  }

  if(family=="gaussian"){
  marg.tau<-inla.tmarginal(function(x){1/x},m.blocknngp$marginals.hyperpar$`Precision for the Gaussian observations`)
  tau.summary<-as.matrix(inla.zmarginal(marg.tau,silent = T))
  marg.sig<-inla.tmarginal(function(x){exp(-x)},m.blocknngp$marginals.hyperpar$`Theta1 for idx`)
  sigma.summary<-as.matrix(inla.zmarginal(marg.sig,silent = T))
  marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.blocknngp$marginals.hyperpar$`Theta2 for idx`)
  phi.summary<-as.matrix(inla.zmarginal(marg.phis,silent = T))

   return(list(res=m.blocknngp,tau.summary=tau.summary,sigma.summary=sigma.summary,phi.summary=phi.summary))
  }else{
    marg.sig<-inla.tmarginal(function(x){exp(-x)},m.blocknngp$marginals.hyperpar$`Theta1 for idx`)
    sigma.summary<-as.matrix(inla.zmarginal(marg.sig,silent = T))
    marg.phis<-inla.tmarginal(function(x){(a+b*exp(x))/(1+exp(x))},m.blocknngp$marginals.hyperpar$`Theta2 for idx`)
    phi.summary<-as.matrix(inla.zmarginal(marg.phis,silent = T))
    return(list(res=m.blocknngp,sigma.summary=sigma.summary,phi.summary=phi.summary))
  }

}



#' Predict spatial effect at new location
#'
#' Draw posterior predictive samples for spatial effect at new locations using fitted NNGP or block-NNGP model. Only new locations are allowed,
#' any location existed in data may cause error.
#' @param n.sample Number of posterior samples needed.
#' @param new.coord A L*2 matrix of new coordinates whose spatial effect needs to be predicted.
#' @param preprocess.data Preprocessed data from \code{\link{preprocess.nngp}} or \code{\link{preprocess.blocknngp}}.
#' @param res Result from \code{\link{nngp.inla}} or \code{\link{blocknngp.inla}}.
#' @param model "nngp" or "blocknngp"; default is "nngp", anything other than "nngp' would be considered as "blocknngp".
#' @param phi.prior 2*1 vector c(a,b) indicates the uniform prior U(a, b) for spatial decay parameter phi.
#' @param mc.cores The number of CPU cores should be used. For Windows user, only 1 core is allowed.
#' @return Function that returns posterior samples.
#' \itemize{
#' \item res - n.sample*L matrix that i_th column contains n.sample posterior samples of i_th location in new.coord.
#' }
#' @export
get.new.spatial<-function(n.sample=1,new.coord ,result, preprocess.data, model="nngp", phi.prior=c(1,30),
                         mc.cores=1 ){
  #### predict new location
  cat("Model:",model)
  require(INLA)
  require(parallel)
  require(pdist)
  if(model=="nngp"){
  M=preprocess.data$M
  train.coord<-preprocess.data$unique.coord
  a<-phi.prior[1]; b<-phi.prior[2]
  pre_Z.inla<-matrix(NA, nr=n.sample,nc=dim(new.coord))

  inla.pos<-inla.posterior.sample(n=n.sample,result$res,num.threads = mc.cores,intern = FALSE,selection = list(idx=0))

    pre.spatial<- function(j="iterations",i="ith location"){
      tt<-as.matrix(pdist(new.coord[i,],train.coord[,1:2]))
      pp<-sort(tt,decreasing = FALSE,index.return=TRUE)
      dist.NN<-pp$x[1:M] #dist to nearest neighbors
      ind.NN<-pp$ix[1:M] #ind of nearest neighbors

      theta1<-inla.pos[[j]]$hyperpar[2]
      theta2<-inla.pos[[j]]$hyperpar[3]

      sigma2<-exp(-theta1)
      phis<-(a+b*exp(theta2))/(1+exp(theta2))

      C_sn<-sigma2*exp(-as.matrix(dist.NN)*phis)
      D<-as.matrix(dist(train.coord[ind.NN,1:2]))
      C_N<-sigma2*exp(-D*phis)
      B.new<-t(C_sn)%*%solve(C_N)
      F.new<-sqrt(sigma2-B.new%*%C_sn)

      if(sigma2-B.new%*%C_sn<0){ print("error:sigma2-B.new%*%C_sn<0")}

      prediction.z<-rnorm(1,B.new%*%inla.pos[[j]]$latent[ind.NN],F.new)
    }

    for (i in 1:dim(new.coord)[1]) {
      cat(paste("\n predict location",i))
      res1<-mclapply(1:n.sample,pre.spatial,i=i,mc.cores = mc.cores)
      pre_Z.inla[,i]<-unlist(res1)
    }
    res<-pre_Z.inla

  }else{

    nb<-preprocess.data$cut.block
    train.coord<-preprocess.data$unique.coord
    neighbor.block<-preprocess.data$neighbor.block
    block.center<-preprocess.data$block.center
    D_bN<-list()
    D_NN<-list()
    D_b<-list()
    for (i in 1:(nb[1]*nb[2])) {
      if(i==1){
        D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, 1:2]))

      }else{
        nn.block<-subset(train.coord,block.ind %in% neighbor.block$NN_ind[i-1,])
        D_bN[[i]]<-as.matrix(pdist(train.coord[train.coord$block.ind==i, 1:2],
                                   nn.block[, 1:2]))
        D_NN[[i]]<-as.matrix(dist(nn.block[,1:2]))
        D_b[[i]]<-as.matrix(dist(train.coord[train.coord$block.ind==i, 1:2]))
      }
    }


    block.nngp.pos<-inla.posterior.sample(n=n.sample,result$res,num.threads = mc.cores,intern = FALSE,selection = list(idx=0))

    block.center<-block.center[neighbor.block$ord,]
    block.center[,1]<-1:(nb[1]*nb[2])

    a<-phi.prior[1]; b<-phi.prior[2]

    predict.spatial<- function(j="iterations",i="ith location"){

      theta1<-block.nngp.pos[[j]]$hyperpar[2]
      theta2<-block.nngp.pos[[j]]$hyperpar[3]

      sigma2<-exp(-theta1)
      phis<-(a+b*exp(theta2))/(1+exp(theta2))

      b.i<-which.min(as.matrix(pdist(new.coord[i,],block.center[,2:3])))
      C_iN<-sigma2*exp(-phis*as.matrix(pdist(new.coord[i,],train.coord[train.coord$block.ind==b.i,1:2])))
      C_N.inv<-solve(sigma2*exp(-phis*D_b[[b.i]]))
      B.i<-C_iN %*% C_N.inv
      F.i<-sigma2-B.i%*%t(C_iN)
      z.i<-rnorm(1, B.i%*% block.nngp.pos[[j]]$latent[which(train.coord$block.ind==b.i)],sqrt(F.i))
      return(z.i)
    }

    pre_Z.block.nngp<-matrix(NA, nr=n.sample,nc=dim(new.coord)[1])
    for (i in 1:dim(new.coord)[1]) {
      cat(paste("\n predict location",i))
      res1<-mclapply(1:n.sample,predict.spatial,i=i,mc.cores = mc.cores)
      pre_Z.block.nngp[,i]<-unlist(res1)
    }
    res<-pre_Z.block.nngp

  }

  return(res)

}


#' Spatial plot
#'
#' Create spatial interpolation plot. Specify result and preprocess.data to get estimated spatial surface of fitted data.
#' @param result Result from \code{\link{nngp.inla}} or \code{\link{blocknngp.inla}}.
#' @param preprocess.data Preprocessed data from \code{\link{preprocess.nngp}} or \code{\link{preprocess.blocknngp}}.
#' @param coord K*2 coordinates matrix. Specify for your customized coordinates, do not need to specify when you want estimated spatial surface of fitted data.
#' @param values 1*k vector of values at corresponding coordinates. Specify for your customized values, do not need to specify when you want estimated spatial surface of fitted data.
#' @param nx Dimension of output grid in x direction.
#' @param ny Dimension of output grid in x direction.
#' @param zlim Adjust the ranges of values.
#' @return Function that returns a figure.
#' @export
spatial.plot<-function( result, preprocess.data,coord=preprocess.data$unique.coord, values=result$res$summary.random$idx$mean,
                       nx=200, ny=200,zlim=NULL){

  require(akima)
  require(fields)

  hh<-coord[,1]
  vv<-coord[,2]
  zz<-values
  sp<-interp(hh,vv,zz,nx=nx,ny=ny)
  image.plot(sp,zlim=zlim)
}

#' Posterior sampling
#'
#' This function generate samples, and functions of those, from an approximated posterior of a fitted model (an inla-object).
#' Please refer to \code{\link[INLA]{inla.posterior.sample}}. The latent spatial effect has name "idx".
#' @import INLA
#' @export
posterior.sample<- inla.posterior.sample
