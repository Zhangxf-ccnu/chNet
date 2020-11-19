
#' @title Differential network analysis by simultaneously considering  changes in gene interactions and gene expression
#' @description The complete procedure for estimating differential network
#'  using chNet. For details, refer to method (Section 2.4 in the main text).
#' @importFrom  parallel mclapply
#' @importFrom  stats pnorm qnorm
#' @importFrom  igraph  graph_from_adjacency_matrix
#' @import Matrix mvtnorm glmnet
#' @usage chNet (X,group,lambar,parallel,nCpus)
#' @param X Data matrix for which the rows represent the samples and the columns represent the genes.
#' @param group Vector which defines two groups of samples under comparison.
#' @param lambar The tuning (threshold) parameter controls the sparsity level.
#' @param parallel Logical value to indicate if the process should be run parallelly in multiple threads, default to FALSE.
#' @param nCpus number of (maximum) cores to use for parallel execution, default to 4.
#' @details This function is implemented to infer diffferential network that satisfy hierarchichal constraints.
#' We first define the differential network as the difference of partial correlations,
#' and develop a new  test statistic to quantify the change of partial correlations.
#' Then the Student's t-test statistic is used to quantify the changes in expression levels of individual genes.
#' Finally, an optimization model is developed to combine the two different types of test statistics so that the
#' estimated differential networks exhibit the hierarchical structures.
#' A closed-formed solution is derived to solve the optimization model.
#' @return
#' \item{\code{diff.edge}}{the adjacency matrix of the estimated differential network.}
#' \item{\code{diff.gene}}{the adjacency matrix of the estimated differentially expressed genes.}
#' \item{\code{Diff.net}}{the estimated differential network over all genes.}
#' \item{\code{Diff.net.connected}}{the estimated differential network over only the connected genes.}
#' @references  Jia-Juan Tu, Le Ou-Yang, Hong Yan, Hong Qin and Xiao-Fei Zhang(2020), Differential network analysis
#'  by simultaneously considering  changes in  gene interactions and gene expression.
#' @author Jia-Juan Tu
#' @seealso { \code{\link{generate.data}}, \code{\link{TCGA.BRCA}}, \code{\link{GSE13159.AML}}}
#' @export
#' @examples # Simulation data
#' data.x= generate.data(p = 100, n = 100, umin = 0.5, umax = 1)
#' result = chNet(data.x$X,data.x$group, lambar = 4, parallel = FALSE, nCpus = 4)
#'
#' # TCGA breast cancer data
#' data("TCGA.BRCA")
#' result = chNet(TCGA.BRCA$X,TCGA.BRCA$group, lambar = 2.85, parallel = FALSE, nCpus = 4)
#' # GSE13159 AML
#' data("GSE13159.AML")
#' result = chNet(TCGA.BRCA$X,TCGA.BRCA$group, lambar = 2.7, parallel = FALSE, nCpus = 4)


chNet <- function(X, group, lambar, parallel = FALSE, nCpus = 4){

  cc = compute.statis(X,group, parallel, nCpus)

  lams.est  <- get.knots(cc$tmu, cc$trho)
  lams.edge.temp =lams.est$int.knots
  diag(lams.edge.temp)=0
  lams.edge  <- (lams.edge.temp >= lambar)*1
  lams.diffgene  <- (lams.est$main.knots >= lambar)*1

  Diff.edge.temp = upper.tri(lams.edge)*lams.edge
  Diff.edge=Diff.edge.temp+t(Diff.edge.temp)
  row.names(Diff.edge) = colnames(X)
  colnames(Diff.edge) = colnames(X)
  Degree.Diff = apply(Diff.edge, 1, sum)
  Diff.graph.connected = Diff.edge[Degree.Diff>0, Degree.Diff>0]
  Diff.net.connected  = graph_from_adjacency_matrix(Diff.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  Diff.net = graph_from_adjacency_matrix(Diff.edge, mode = "undirected", weighted = TRUE, diag = FALSE)

  return(list(diff.edge=lams.edge, diff.gene=lams.diffgene,
              Diff.net=Diff.net,Diff.net.connected=Diff.net.connected))
}



lasso.estimator <- function(X, lam1,  parallel = FALSE, nCpus = 4){


  n = dim(X)[1]
  p = dim(X)[2]

  centered.X = scale(X, scale = FALSE)
  Sigma = cov(centered.X)
  DSigma = diag(Sigma)

  lam2 = sqrt(DSigma*log(p)/n)
  beta = array(0, dim=c(p,p,length(lam1)))
  res = array(0,dim=c(n,p,length(lam1)))

  wrapper = function(i){
    fit=glmnet(centered.X[,-i], centered.X[,i], lambda= lam1*lam2[i])
    fit$beta=as.matrix(fit$beta)
    if(ncol(fit$beta)<length(lam1)){
      tmp = matrix(0,nrow = nrow(fit$beta),ncol = length(lam1))
      tmp[,1:ncol(fit$beta)]=fit$beta
      tmp[,ncol(fit$beta):length(lam1)] = fit$beta[,ncol(fit$beta)]
      fit$beta = tmp
    }

    if(i==1){
      beta[2:p,i,]=fit$beta
    }else if(i==p){
      beta[1:(p-1),i,]=fit$beta
    }else{
      beta[1:(i-1),i,]=fit$beta[1:(i-1),]
      beta[(i+1):p,i,]=fit$beta[i:nrow(fit$beta),]
    }
    res[,i,] = matrix(rep(centered.X[,i],length(lam1)),ncol = length(lam1)) - centered.X[,-i]%*%fit$beta
    out = list(beta = beta[,i,], res = res[,i,])
    return(out)
  }


  if(parallel){
    fit =mclapply(1:p, wrapper, mc.cores=nCpus)
  }else{
    fit =lapply(1:p, wrapper)
  }


  for(i in 1:p){
    beta[,i,]=fit[[i]]$beta
    res[,i,]=fit[[i]]$res
  }


  r.tilde = array(0, dim=c(p,p,length(lam1)))
  r.hat = array(0, dim=c(p,p,length(lam1)))

  for(k in seq(length(lam1))){
    r.tilde[,,k] = cov(res[,,k])*(n-1)/(n)
    r.hat[,,k] =r.tilde[,,k] + diag(diag(r.tilde[,,k]))%*%beta[,,k] + t(beta[,,k])%*%diag(diag(r.tilde[,,k]))
  }

  out = list(beta = beta, res = res, r.tilde = r.tilde, r.hat = r.hat)
  return(out)
}

## calculate test statistics
compute.statis<- function(X, group, parallel = FALSE, nCpus = 4 ) {

  n <- nrow(X)
  p <- ncol(X)

  lam1 = seq(40,4)/20
  nlam1=length(lam1)

  stopifnot(length(group) == n)

  stopifnot(table(group) > 3) # need at least 2 observations in each class

  X1 = X[group==unique(group)[1],]
  X2 = X[group==unique(group)[2],]

  n1 <- nrow(X1)
  n2 <- nrow(X2)

  #step 1  calculate the t-statistical of mu
  xbar1 <- colMeans( X1)
  xbar2 <- colMeans( X2)
  s1 <- apply( X1, 2, sd)
  s2 <- apply( X2, 2, sd)

  tmu <- (xbar1-xbar2) / sqrt(s1^2/n1+ s2^2/n2)


  #step 2  calculate the t-statistical of rho

  fit1 = lasso.estimator(X1, lam1, parallel, nCpus)
  fit2 = lasso.estimator(X2, lam1, parallel, nCpus)

  WW = array(0, dim=c(p,p,nlam1))
  score = rep(0,nlam1)
  for(k in seq(nlam1)){
    # print(k)
    Dr1 = diag((diag(fit1$r.tilde[,,k]))^(-0.5))
    Dr2 = diag((diag(fit2$r.tilde[,,k]))^(-0.5))
    T1 = Dr1%*%fit1$r.hat[,,k]%*%Dr1
    T2 = Dr2%*%fit2$r.hat[,,k]%*%Dr2
    rho1 = T1*((abs(T1)>=2*sqrt(log(p)/n1))*1)
    rho2 = T2*((abs(T2)>=2*sqrt(log(p)/n2))*1)
    WW[,,k] = (T1-T2)/sqrt((((1-rho1^2)^2)/n1)+(((1-rho2^2)^2)/n2))
    diag(WW[,,k] ) = 0
    WW[,,k] = WW[,,k]*upper.tri(WW[,,k])
    for (l in seq(1,10)){
      score[k] = score[k] + (sum(abs(WW[,,k])>qnorm(1-l*(1-pnorm(sqrt(log(p))))/10))/(l*(p^2-p)*(1-pnorm(sqrt(log(p))))/10)-1)^2
    }
  }

  kmin = which(score==min(score))

  if(length(kmin)>1){
    k.min=kmin[1]
  }else{k.min=kmin}

  W.opt = WW[,,k.min]

  W.opt =as.matrix(  W.opt )
  diag(W.opt) = 0

  return(list(tmu=tmu, trho=W.opt))

}

# A closed-formed solution is derived to solve the optimization model.
get.knots <- function(w, z) {

  p <- length(w)
  #cat("p=",p)
  #cat("nrow(z) =",nrow(z) )
  if (nrow(z) != p) browser()
  stopifnot(nrow(z) == p, ncol(z) == p)
  aw <- abs(w)
  diag(z) <- rep(0, p)
  az <- abs(z)

  main.knots <- pmax(aw, (aw + apply(az, 1, max)) / 2)
  int.knots <- matrix(NA, p, p)
  for (j in seq(p)) {
    for (k in seq(p)) {
      if (j == k) next
      int.knots[j, k] <- max(aw[j] - sum(pmax(az[j, -j] - az[j, k], 0)), 0) / 2
    }
  }
  int.knots <- pmin(az, az / 2 + int.knots)
  int.knots <- pmax(int.knots, t(int.knots)) # use max(lam_jk, lam_kj)

  return(list(main.knots=main.knots, int.knots=int.knots))

}

