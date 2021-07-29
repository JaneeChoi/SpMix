#' @importFrom LogConcDEAD mlelcd
#' @importFrom mvtnorm dmvnorm
#'
#'
#' @title Estimates a Semiparametric Mixture Density for multi-dimensional data
#'
#' @description \code{sp.mix.multi} returns a semiparametric mixture density estimates for given multi-dimensional z, which are the probit-transformed p-values.
#'
#' @param z Matrix which each column indicates z, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm.
#' Optimization stops if maximum absolute difference of current and previous gamma value is smaller than tol. (default: 5.0e-6)
#' @param max.iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param mono If TRUE, monotone. (default: TRUE)
#'
#' @return Parametrization of f(x) in terms of hyperplanes and function
#'   evaluations y = log(f(x)) \item{aOpt, bOpt}{Analytically normalized
#'   parameters of f(x).} \item{logLike}{Log-likelihood of f(x)} \item{y}{Vector
#'   with values y_i = log(f(X_)) of the normalized density (\eqn{logLike =
#'   \sum(y_i)}).} \item{aOptSparse, bOptSparse}{Sparse parametrization
#'   normalized on the integration grid.}
#'
#' @export
sp.mix.multi <- function(z, tol = 5e-6, max.iter = 30, mono = TRUE)
  # FOR MULTIVARIATE CASE ONLY
{

  z <- as.matrix(z)
  n <- dim(z)[1]

  ## Initial step: to fit normal mixture
  nmEM <- normal.mixture(z)
  p.0 <- nmEM$p.0
  mu.0 <- nmEM$mu.0
  sig.0 <- nmEM$Sigma.0
  f1.tilde <- dmvnorm(z, nmEM$mu.1, nmEM$Sigma.1)
  gam <- f <- rep(0, n)

  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3)|((k < max.iter) & (!converged)) ) {
    k <- k + 1
    ## E-step
    tmp <- p.0*dmvnorm(z, mu.0, sig.0)
    new.f <- tmp + (1-p.0)*f1.tilde
    new.gam <- tmp/new.f
    if(mono) new.gam <- MonotoneFDR(z, new.gam)

    ## M-step
    sum.gam <- sum(new.gam)
    new.mu.0 <- as.vector(t(z)%*%new.gam)/sum.gam
    dev <- t(t(z)-new.mu.0)*sqrt(new.gam)
    new.sig.0 <- t(dev)%*%dev/sum.gam
    new.p.0 <- mean(new.gam)
    new.f.0 <- dmvnorm(z, new.mu.0, new.sig.0)
    weight <- 1 - new.gam
    new.f1.tilde <- rep(0, n)
    which.z <- (new.gam <= .9)
    lcd <- LogConcDEAD::mlelcd(z[which.z,], w = weight[which.z]/sum(weight[which.z]))
    new.f1.tilde[which.z] <- exp(lcd$logMLE)

    ## Update
    which.gam <- (new.gam <= 0.9)*(new.gam >= 0.01)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in mdfdr fit = ", round(diff, 5), "\n")
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde
    gam <- new.gam
    f <- new.f
  }

  res <- list(p.0 = p.0, mu.0 = mu.0, tau.0 = sig.0,
              f1.hat = f1.tilde, f = f, localfdr = gam, iter = k)

  return(res)
}

normal.mixture <- function(z, tol = 5e-3, max.iter = 10)
{

  k <- 0; diff <- 100
  z <- as.matrix(z)
  if ( dim(z)[2] == 1 ) m.dist <- z/sd(z) else m.dist <- mahalanobis(z, center = rep(0, dim(z)[2]), cov = cov(z))

  p.0 <- mean(m.dist <= 1.65)
  mu.0 <- rep(0, dim(z)[2])
  sig.0 <- diag(1, dim(z)[2])
  f.0 <- dmvnorm(z, mean = mu.0, sigma = sig.0)

  mu.1 <- apply(z[m.dist > 1.65,], 2, mean)
  sig.1 <- cov(z[m.dist > 1.65,])
  f.1 <- dmvnorm(z, mean = mu.1, sigma = sig.1)

  while ( (k < 3)|((k < max.iter) & (diff > tol)) ) {
    k <- k + 1

    ## E-step
    term1 <- p.0*f.0
    term2 <- term1 + (1-p.0)*f.1
    gam <- term1/term2

    ## M-step
    new.p.0 <- mean(gam)
    new.mu.0 <- as.vector(t(z)%*%gam)/sum(gam)
    dev <- (z - new.mu.0)*sqrt(gam)
    new.sig.0 <- t(dev)%*%dev/sum(gam)
    f.0 <- dmvnorm(z, mean = new.mu.0, sigma = new.sig.0)
    new.mu.1 <- as.vector(t(z)%*%(1 - gam))/sum(1 - gam)
    dev <- (z - new.mu.1)*sqrt(1 - gam)
    new.sig.1 <- t(dev)%*%dev/sum(1 - gam)
    f.1 <- dmvnorm(z, mean = new.mu.1, sigma = new.sig.1)

    ## Update
    diff <- max(abs((new.mu.0 - mu.0)),
                abs((new.sig.0 - sig.0)),
                abs((new.mu.1 - mu.1)),
                abs((new.sig.1 - sig.1)),
                abs((new.p.0 - p.0)))
    p.0 <- new.p.0
    mu.0 <- new.mu.0
    sig.0 <- new.sig.0
    mu.1 <- new.mu.1
    sig.1 <- new.sig.1
  }
  return(list(iter = k, p.0 = p.0, mu.0 = mu.0, Sigma.0 = sig.0, mu.1 = mu.1, Sigma.1 = sig.1))
}



NE <- function(x, X)
{
  n <- nrow(X)
  p <- ncol(X)
  xx <- matrix(x, nrow = n, ncol = p, byrow = TRUE)
  ne.ind <- apply(1*(X >= xx), 1, prod)

  return((1:n)[ne.ind == 1])
}


MonotoneFDR <- function(z, fdr)
{
  n <- nrow(z)
  MFDR <- numeric(n)
  for (i in 1:n) {
    MFDR[i] <- max(fdr[NE(z[i,], z)])
  }

  return(MFDR)
}
