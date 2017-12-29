#' Applying Lauritzen's formula on covariance matrix
#'
#' This function applies the lauritzen's formula to a covariance matrix
#' to take into account a decomposable graph structure. It uses the provided
#' covariance matrix and the provided graph to compute the covariance matrix
#' that perfectly fits the set of conditional independencies
#' encoded by the graph's cliques and seperators.
#'
#' @param A An initial covariance Matrix.
#' @param C A list of numeric vectors which are the graph's maximal cliques.
#' @param S A list of numeric vectors which are the graph's minimal seperators.
#' @param nS A numeric vector which indicates the multiplicity of each
#' minimal seperator. It must have the same length as S.
#'
#' @return This function returns a covarianece matrix which takes into
#' account the provided graph structure and the initial covariance matrix.
#'
#' @examples
#' A <- matrix(0,5,5)
#' diag(A) <- 1:5
#' diag(A[-1,-5]) <- 1:4
#' diag(A[-5,-1]) <- 1:4
#' cliques <- list(c(1,2),c(2,3), c(3,4),c(4,5))
#' separators <- list(c(2),c(3),c(4))
#' nS <- c(1,1,1)
#' Anew <- graphSigma(A, cliques, separators, nS)
#'
#' @export
graphSigma <- function(A, C, S, nS){
  d <- nrow(A)
  Q <- matrix(0,d,d)

  for(clique in C){
    Ac <- A[clique,clique]
    Q[clique,clique] <- Q[clique,clique] + solve(Ac)
  }
  k=0
  for(sep in S){
    k <- k+1
    As <- A[sep,sep]
    Q[sep,sep] <- Q[sep,sep] - solve(As)*nS[k]
  }
  res <- solve(Q)
  res
}


#' Measure the association degree between a graph and a covariance Matrix
#'
#' This function computes the association degree between a covariabce matrix
#' and a graph. The computed metric relies only on the correspondance between
#' the zeros in the inverted covariance matix and the set of conditional independencies.
#'
#' @param A A covariance Matrix.
#' @param C A list of numeric vectors which specify the graph's maximal cliques
#'
#' @return This function returns a numeric value which measures the association degree between the covariance matrix and the graph.
#'
#' @examples
#' A <- matrix(0,5,5)
#' diag(A) <- 1:5
#' diag(A[-1,-5]) <- 1:4
#' diag(A[-5,-1]) <- 1:4
#' cliques <- list(c(1,2),c(2,3), c(3,4),c(4,5))
#' d1 <- graphMatrixAssoc(A,C)
#'
#' @export
graphMatrixAssoc <- function(A, C){
  Q <- solve(A)
  d <- nrow(A)

  nadjmat = matrix(1,d,d)
  for(clique in C){
    nadjmat[clique,clique] <- 0
  }
  nQ <- sqrt(sum(diag(Q%*%Q)))
  dist <- sum(abs(Q)*nadjmat)/nQ
  dist
}

#' The posterior probability
#'
#' This function calculates the posterior probability that each observation belongs to each of the mixture components.
#'
#' @param x A matrix of observations with d rows (number of variables)
#' and n columns (number of observations).
#' @param p A numeric vector of proportions (one for each mixture component).
#' @param mu A list of numeric vectors (which are the mean of each mixture component).
#' @param sigma A list of covariance matrices (one for each mixture component).
#'
#' @return This function returns a matrix that contains the posterior probabilities
#' @importFrom mvtnorm dmvnorm
#'
#' @export

computeTau <- function(x, p, mu, sigma){
  n <- ncol(x)
  d <- nrow(x)
  K <- length(p)
  Tho <- matrix(0,n,K)
  for (i in 1:n){
    s =0;
    for(j in 1:K){
      Tho[i,j]=p[j]*dmvnorm(x[,i],mu[[j]],sigma[[j]])
      s = s+Tho[i,j]
    }

    if(s!=0){
      Tho[i,]=Tho[i,]/s
    }
  }
  Tho
}

#' Estimates the Gaussian mixture parameters using GEM algorithm
#'
#' This function estimates the Gaussian mixture parameters using the GEM
#' algorithm. The mixture components number is supposed to be known.
#' This function uses an initial parameters guess and a set of associated
#' graphs to iteratively estimate the parameters.
#'
#' @param x A matrix of observations with d rows (number of variables) and n columns (number of observations).
#' @param Nit An integer that specifies the iterations number. If a negative value is
#' specified then a stopping rule is used.
#' @param p0 A numeric vector of initial proportions.
#' @param mu0 A list of numeric vectors which are the initial mean of each
#' mixture component.
#' @param sigma0 A list of covariance matrices (Algorithm initialization)
#' @param G A list of K graphs associated with the mixture components. Each
#' graph is a named list with three items : Cliques list, Sperators list and
#' Sperators multiplicities. If an empty list is provided then the classical
#' EM algorithm is used.
#'
#' @return This function returns a named list containing the following results:
#' 'p': A vector of estimated proportions.
#' 'mu': A list of estimated mean vectors.
#' 'sigma' : A list of estimated covarianece matrices which takes into
#' account the provided graph structures.
#' 'tau' : An n by K matrix that contains the posterior probabilities that each observation belongs to each mixture component.
#' 'NbIterations' ; The number of performed iterations.
#'
#' @export
gemEstimator <- function(x, Nit, p0, mu0, sigma0, G=list()){
  n <- ncol(x)
  d <- nrow(x)
  K <- length(p0)
  Tau <- computeTau(x,p0,mu0,sigma0);

  mu <- mu0;sigvect <- sigma0;p <- p0;

  applygem <- length(G)>0
  useStop <- FALSE

  if(Nit >0)
    Nitmax <- Nit
  else{
    Nitmax <- 100
    useStop <- TRUE
  }
  NbIterations = 0
  epsStop = 1e-4
  notok = TRUE
  while(NbIterations < Nitmax & notok){
    d0 <- normTheta(p,mu,sigvect)
    previousSig = sigvect
    previousMu = mu
    previousp = p
    for( k in 1:K){
      nk <- sum(Tau[,k])
      p[k] <- nk/n
      mu[[k]] <- sapply(1:d, function(idx){sum(Tau[,k]*t(x[idx,]))/nk})
      S<- matrix(0,d,d)
      for(i in 1:n){
        S <- S + Tau[i,k]*((x[,i]-mu[[k]])%*%(t(x[,i])-mu[[k]]))
      }
      if(applygem)
        sigvect[[k]] <- graphSigma(S/nk, C = G[[k]][["Cliques"]],S = G[[k]][["Separators"]], nS = G[[k]][["multiplicities"]])
      else
        sigvect[[k]] <- S/nk
    }
    Tau <- computeTau(x,p,mu, sigvect)
    NbIterations <- NbIterations +1
    d1 <- normTheta(p - previousp, mapply('-',mu, previousMu,SIMPLIFY = F),
                    mapply('-',sigvect, previousSig,SIMPLIFY = F) )
    notok <- (d1 > epsStop*(1+d0)) & useStop
  }
  res <- list(
    p = p,
    mu = mu,
    sigma = sigvect,
    Tau = Tau,
    NbIterations = NbIterations
  )
  res
}

#' Computes the norm of a set of parameters.
#'
#' This function computes the norm of the Gaussian mixture parameters set as the
#' sum of the respective norms of each element.
#'
#' @param p A numeric vector of size K.
#' @param mu A list of K numeric vectors of size d.
#' @param sigma A list of K square matrices of size dxd
#'
#' @return This function returns a numeric value which is the norm of the provided
#' parameters.


normTheta <- function(p, mu, sigma){
  d <- sum(abs(p)) + sum(sapply(mu, norm, type="2")) + sum(sapply(sigma, function(sig){sqrt(sum(diag(sig%*%sig)))}))
  d
}
