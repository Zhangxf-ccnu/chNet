
#' @title Generate simulated data
#' @description The complete procedure for generating simulated data. For details, refer to
#' simulation study (Supplementary Section 3.3.1 ).
#' @importFrom igraph as_adjacency_matrix sample_pa
#' @importFrom  MASS mvrnorm
#' @importFrom  stats rbinom rnorm runif sd cov
#' @usage generate.data (p, n, umin, umax)
#' @param p The number of genes.
#' @param n The sample size.
#' @param umin  The lower limits of the edge values.
#' @param umax  The upper limits of the edge values.
#' @details The function is used to generate the gene expression datasets.
#' @return
#' \item{\code{X}}{ A matrix of sample  matrice  (\eqn{2n \times p}) from two different conditions.}
#' \item{\code{group}}{ A matrix of sample label matrice  (\eqn{2n \times 1}) from two different conditions.}
#' \item{\code{rho}}{A list (length = 2) of the partical coefficients matrices (\eqn{p \times p}).}
#' \item{\code{hubgene}}{A set hub genes.}
#'
#'
#' @references  Jia-Juan Tu, Le Ou-Yang, Hong Yan, Hong Qin and Xiao-Fei Zhang (2021), Differential network analysis
#'  by simultaneously considering  changes in  gene interactions and gene expression.
#' @author Jia-Juan Tu
#' @seealso { \code{\link{chNet}}, \code{\link{TCGA.BRCA}}, \code{\link{GSE13159.AML}}}
#' @export
#' @examples
#' # Simulation data
#' data.x = generate.data(p = 100, n = 100, umin = 0.5, umax = 1)
generate.data <-  function(p = 100, n = 100, umin = 0.5, umax = 1){

  A = as_adjacency_matrix(sample_pa(p, directed = FALSE), type = "both", sparse=FALSE)


  alpha=0.1
  W = matrix(runif(p*p, min = umin, max = umax)*(2*rbinom(p*p, 1, 0.5) - 1), p, p)
  W1 = A*W
  W1 = W1 + t(W1)
  W2 = W1

  #Step 1: select the hub genes
  diffnum = floor(alpha*p)
  hub_ind_diff = sample(1:diffnum, 3)

  hub_ind_nondiff = sample((diffnum+1):p, 1)
  hubgene = c(hub_ind_diff, hub_ind_nondiff )

  #########Step 2: create the differential network

  ######### hub genes are differentially expressed genes
  for(i in 1: 3){
    Indx1 = sample(1:p, diffnum )
    ind_value = runif(length(Indx1), min = umin, max = umax)*(2*rbinom(length(Indx1), 1, 0.5) - 1)

    indy1 = sample(1:2,1)

    if (indy1==1){
      W1[hub_ind_diff[i],Indx1] = ind_value
      W1[Indx1,hub_ind_diff[i]] = ind_value
    }else{
      W2[hub_ind_diff[i],Indx1] = ind_value
      W2[Indx1,hub_ind_diff[i]] = ind_value
    }
  }

  ######### hub gene is a non-differentially expressed genes

  Indx2 = sample(1:diffnum, diffnum)
  ind_value = runif(length(Indx2), min = umin, max = umax)*(2*rbinom(length(Indx2), 1, 0.5) - 1)

  indy2 = sample(1:2,1)

  if (indy2==1){
    W1[hub_ind_nondiff,Indx2] = ind_value
    W1[Indx2,hub_ind_nondiff] = ind_value
  }else{
    W2[hub_ind_nondiff,Indx2] = ind_value
    W2[Indx2,hub_ind_nondiff] = ind_value

  }


  eigen.min = min(eigen(W1)$values,eigen(W2)$values)

  W1 = W1 + (abs(eigen.min) + 0.1)*diag(p)
  W2 = W2 + (abs(eigen.min) + 0.1)*diag(p)


  D1 = diag((diag(W1))^(-0.5))
  R1 = D1%*%W1%*%D1
  D2 = diag((diag(W2))^(-0.5))
  R2 = D2%*%W2%*%D2


  diag1 = rep(1,p)
  diag2 = rep(1,p)
  diag2[1:floor(alpha*p)] = 4

  Omega1 =  diag(sqrt(diag1))%*%R1%*%diag(sqrt(diag1))
  Omega2 =  diag(sqrt(diag2))%*%R2%*%diag(sqrt(diag2))

  mu1=rep(0,p)
  mu2=mu1
  mu2[1:floor(alpha*p)]=2
  X1 = mvrnorm(n, mu1, solve(Omega1))
  X2 = mvrnorm(n, mu2, solve(Omega2))

  X = rbind(X1, X2)
  group = c(rep(1,n), rep(2,n))

  Omega = list()
  Omega[[1]] = Omega1
  Omega[[2]] = Omega2


  rho = list()
  rho[[1]] = - diag((diag(Omega1))^(-1/2))%*%Omega1%*%diag((diag(Omega1))^(-1/2))
  rho[[2]] = - diag((diag(Omega2))^(-1/2))%*%Omega2%*%diag((diag(Omega2))^(-1/2))

  result = list(X = X, group = group, rho = rho, hubgene = hubgene)

  result
}
