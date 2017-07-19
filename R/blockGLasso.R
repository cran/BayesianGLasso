#' Block Gibbs sampler function
#'
#' Blockwise sampling from the conditional distribution of a permuted column/row
#' for simulating the posterior distribution for the concentration matrix specifying
#' a Gaussian Graphical Model
#' @param X Data matrix
#' @param iterations Length of Markov chain after burn-in
#' @param burnIn Number of burn-in iterations
#' @param lambdaPriora Shrinkage hyperparameter (lambda) gamma distribution shape
#' @param lambdaPriorb Shrinkage hyperparameter (lambda) gamma distribution scale
#' @param verbose logical; if TRUE return MCMC progress
#' @details Implements the block Gibbs sampler for the Bayesian graphical lasso
#' introduced in Wang (2012). Samples from the conditional distribution of a 
#' permuted column/row for simulating the posterior distribution for the concentration 
#' matrix specifying a Gaussian Graphical Model
#' @return 
#' \item{Sigma}{List of covariance matrices from the Markov chain}
#' \item{Omega}{List of concentration matrices from the Markov chains}
#' \item{Lambda}{Vector of simulated lambda parameters}
#' @author Patrick Trainor (University of Louisville)
#' @author Hao Wang
#' @references Wang, H. (2012). Bayesian graphical lasso models and efficient 
#' posterior computation. \emph{Bayesian Analysis, 7}(4). <doi:10.1214/12-BA729> .
#' @examples
#' \donttest{
#' # Generate true covariance matrix:
#' s<-.9**toeplitz(0:9)
#' # Generate multivariate normal distribution:
#' set.seed(5)
#' x<-MASS::mvrnorm(n=100,mu=rep(0,10),Sigma=s)
#' blockGLasso(X=x)
#' }
#' # Same example with short MCMC chain:
#' s<-.9**toeplitz(0:9)
#' set.seed(6)
#' x<-MASS::mvrnorm(n=100,mu=rep(0,10),Sigma=s)
#' blockGLasso(X=x,iterations=100,burnIn=100)
#' @export
blockGLasso<-function(X,iterations=2000,burnIn=1000,lambdaPriora=1,lambdaPriorb=1/10,
                      verbose=TRUE)
{
  # Total iterations:
  totIter<-iterations+burnIn 
  
  # Sum of product matrix, covariance matrix, n
  S<-t(X)%*%X
  Sigma=stats::cov(X)
  n=nrow(X)
  
  # Concentration matrix and it's dimension:
  Omega<-MASS::ginv(Sigma)
  p<-dim(Omega)[1]
  
  # Indicator matrix and permutation matrix for looping through columns & rows ("blocks")
  indMat<-matrix(1:p**2,ncol=p,nrow=p)
  perms<-matrix(NA,nrow=p-1,ncol=p)
  permInt<-1:p
  for(i in 1:ncol(perms))
  {
    perms[,i]<-permInt[-i]
  }
  
  # Structures for storing each MCMC iteration:
  SigmaMatList<-OmegaMatList<-list()
  lambdas<-rep(NA,totIter)

  # Latent tau:
  tau<-matrix(NA,nrow=p,ncol=p)
  
  # Gamma distirbution posterior parameter a:
  lambdaPosta<-(lambdaPriora+(p*(p+1)/2))
  
  # Main block sampling loop:
  for(iter in 1:totIter)
  {
    # Gamma distirbution posterior parameter b:
    lambdaPostb<-(lambdaPriorb+sum(abs(c(Omega)))/2)
    # Sample lambda:
    lambda<-stats::rgamma(1,shape=lambdaPosta,scale=1/lambdaPostb)
    
    OmegaTemp<-Omega[lower.tri(Omega)]
    #cat("Omega Temp min=",min(OmegaTemp)," max=",max(OmegaTemp),"\n")
    
    # Sample tau:
    rinvgaussFun<-function(x)
    {
      x<-ifelse(x<1e-12,1e-12,x)
      #cat("lambda=",lambda," mu=",x,"\n")
      return(statmod::rinvgauss(n=1,mean=x,shape=lambda**2))
    }
    tau[lower.tri(tau)]<-1/sapply(sqrt(lambda**2/(OmegaTemp**2)),rinvgaussFun)
    tau[upper.tri(tau)]<-t(tau)[upper.tri(t(tau))]
    
    # Sample from conditional distribution by column:
    for(i in 1:p)
    {
      tauI<-tau[perms[,i],i]
      Sigma11<-Sigma[perms[,i],perms[,i]]
      Sigma12<-Sigma[perms[,i],i]
      S21<-S[i,perms[,i]]
      Omega11inv<-Sigma11-Sigma12%*%t(Sigma12)/Sigma[i,i]
      Ci<-(S[i,i]+lambda)*Omega11inv+diag(1/tauI)
      CiChol<-chol(Ci)
      mui<-solve(-Ci,S[perms[,i],i])
      # Sampling:
      beta<-mui+solve(CiChol,stats::rnorm(p-1))
      
      # Replacing omega entries
      Omega[perms[,i],i]<-beta
      Omega[i,perms[,i]]<-beta
      gamm<-stats::rgamma(n=1,shape=n/2+1,rate=(S[1,1]+lambda)/2)
      Omega[i,i]<-gamm+t(beta) %*% Omega11inv %*% beta
      
      # Replacing sigma entries
      OmegaInvTemp<-Omega11inv %*% beta
      Sigma[perms[,i],perms[,i]]<-Omega11inv+(OmegaInvTemp %*% t(OmegaInvTemp))/gamm
      Sigma[perms[,i],i]<-Sigma[i,perms[,i]]<-(-OmegaInvTemp/gamm)
      Sigma[i,i]<-1/gamm
    }
    if(iter %% 100 ==0)
    {
      cat("Total iterations= ",iter, "Iterations since burn in= ", 
          ifelse(iter-burnIn>0,iter-burnIn,0), "\n")
    }
    
    # Save lambda:
    lambdas[iter]<-lambda
    
    # Save Sigma and Omega:
    SigmaMatList[[iter]]<-Sigma
    OmegaMatList[[iter]]<-Omega

  }
  list(Sigmas=SigmaMatList,Omegas=OmegaMatList,lambdas=lambdas)
}

