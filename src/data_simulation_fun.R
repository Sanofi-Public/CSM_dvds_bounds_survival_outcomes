# Functions to generate different datasets


#' Logistic function
#'
#' @param x A number or a vector
#'
#' @return The evaluation of logistic(x).
#' 
#' @keywords internal
logistic <- function(x) {
  return(1 / (1 + exp(-x)))
}


#' Function to generate T0 and T1
#'
#' @param X Matrix of observed confounders
#' @param U Matrix of unobserved confounders
#' @param a Treatment allocation (either 1 for treated, or 0 for control)
#' @param n Sample size
#' @param p.X Dimension of X
#' @param p.U Dimension of U
#'
#' @return A vector of size n
#' 
#' @keywords internal
T.fun <- function(X, U, a, n, p.X, p.U) {
  
  Unif.a <- runif(n=n, min=0, max=1)
  beta.X <- log(runif(n=p.X, min=1.1, max=1.3))
  beta.U <- log(runif(n=p.U, min=0.5, max=0.8))

  ind.nb <- 1
  shape <- 1.8
  fact <- 10
  fact.exp <- 0.95
  const <- 5  # Increase this value to increase the average treatment effect
  lamb <- fact.exp / fact**shape
  scale <- lamb * exp(log(const) * a + beta.X %*% X[ind.nb, ] + beta.U %*% U[ind.nb, ])
  
  T.distrib.fun <- function(x) {return(exp(-lamb * x**shape * exp(log(const) * a + beta.X %*% t(X) + beta.U %*% t(U))))}
  T.gen <- fact * (- log(Unif.a) / (fact.exp * exp(log(const) * a + beta.X %*% t(X) + beta.U %*% t(U))))**(1/shape)
  return(list(T.gen=T.gen, T.distrib.fun=T.distrib.fun))
}


#' Function to generate C0 and C1
#'
#' @param X Matrix of observed confounders
#' @param U Matrix of unobserved confounders
#' @param a Treatment allocation (either 1 for treated, or 0 for control)
#' @param n Sample size
#' @param p.X Dimension of X
#' @param p.U Dimension of U
#' @param shape Shape parameter of a Weibull distribution
#'
#' @return A vector of size n
#' 
#' @keywords internal
C.fun <- function(X, U, a, n, p.X, p.U, shape) {
  
  Unif.a <- runif(n=n, min=0, max=1)
  beta.X <- runif(n=p.X, min=2, max=2.5)
  # Remove the effect of U
  beta.U <- rep(0, p.U)
  
  ind.nb <- 1
  fact <- 10
  fact.exp <- 0.95
  const <- 20  # Increase this value to increase the average treatment effect
  lamb <- fact.exp / fact**shape
  scale <- lamb * exp(log(const) * a + beta.X %*% X[ind.nb, ] + beta.U %*% U[ind.nb, ])
  
  T.distrib.fun <- function(x) {return(exp(-lamb * x**shape * exp(log(const) * a + beta.X %*% t(X) + beta.U %*% t(U))))}
  T.gen <- fact * (- log(Unif.a) / (fact.exp * exp(log(const) * a + beta.X %*% t(X) + beta.U %*% t(U))))**(1/shape)
  return(list(C.gen=T.gen, C.distrib.fun=T.distrib.fun))
}


#' Function to generate the simulated data
#'
#' @param p.X Dimension of X
#' @param p.U Dimension of U
#' @param n Sample size
#' @param alpha Mixture parameter
#' @param Beta Coefficient of X in the definition of U|X=x (matrix of size p.X*p.U)
#' @param delta Parameter for the nominal propensity score on the observational data
#' @param gamma.data True confounding strength
#' @param sigma.G Standard deviation of the Gaussian distribution G in the definition of U|X=x
#' @param shape Shape/scale parameter for the Weibull distribution of C
#' @param inform.cens Whether informative censoring should be used or not
#' @param verbose A boolean indicating whether to print Pearson correlations or not
#'
#' @return A list containing the generated data.
#' @export
#'
#' @examples
genData <- function(p.X=5, p.U=3, n=2000, alpha=0.2,
                    Beta=matrix(runif(p.X*p.U, min=0, max=1),
                                nrow=p.X, ncol=p.U, byrow=TRUE),
                    delta=runif(n=p.X, min=-0.5, max=0.5),
                    gamma.data=5,
                    sigma.G=1,
                    shape=0.354,
                    inform.cens=FALSE,
                    verbose=FALSE) {
  
  if (p.X < 1) {
    stop("p.X should be at least 1")
  }
  
  if (length(delta) != p.X) {
    stop("length(delta) and p.X should be equal")
  }
  
  # Draw X: n vectors of size p.X,
  # uniformly between (-1, 1)
  X <- matrix(runif(n=p.X*n, min=-1, max=1), nrow=n)
  
  # Draw U conditionally on X: n vectors of size p.U
  
  # Mean to generate U|X (p.U*n)
  U.mu <- (1 - alpha) * t(Beta) %*% t(X)
  
  # Be careful, vapply is applied by column!
  U <- t(matrix(vapply(X=t(U.mu), FUN=rnorm, FUN.VALUE=numeric(1),
                       n=1, sd=alpha*sigma.G), nrow=p.U, byrow=TRUE))
  
  e.X <- c(logistic(delta %*% t(X) + 0.5))

  l <- e.X / (e.X + (1 - e.X) * gamma.data)
  u <- e.X / (e.X + (1 - e.X) / gamma.data)
  
  t.X <- qnorm((e.X - l) / (u - l), mean=colMeans(U.mu), sd=alpha*sigma.G/sqrt(p.U))
  
  e.XU <- ifelse(rowMeans(U) > t.X, l, u)
  
  A <- rbinom(n=n, size=1, prob=e.XU)
  
  # Generate T1 and T0
  T1.list <- T.fun(X=X, U=U, a=1, n=n, p.X=p.X, p.U=p.U)
  T0.list <- T.fun(X=X, U=U, a=0, n=n, p.X=p.X, p.U=p.U)
  T1 <- T1.list$T.gen
  T0 <- T0.list$T.gen
  
  T01 <- ifelse(A == 1, T1, T0)
  
  # Generate C (independent censoring, can be changed for conditional independent censoring)
  if (inform.cens) {
    C1 <- C.fun(X=X, U=U, a=1, n=n, p.X=p.X, p.U=p.U, shape=shape)$C.gen
    C0 <- C.fun(X=X, U=U, a=0, n=n, p.X=p.X, p.U=p.U, shape=shape)$C.gen
    C <- ifelse(A == 1, C1, C0)
  } else {
    C <- rweibull(n=n, shape=shape, scale=10)
  }
  
  T.obs <- pmin(T01, C)
  
  c.rate <- mean(C < T01)
  message("Censoring rate = ", c.rate)
  
  return(list(X=X, U=U, T1=T1, T0=T0,
              T1.distrib.fun=T1.list$T.distrib.fun,
              T0.distrib.fun=T0.list$T.distrib.fun,
              T01=T01, C=C, T.obs=T.obs, A=A))
}
