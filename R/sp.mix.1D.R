#' @importFrom logcondens activeSetLogCon
#'
#' @title LocalFDR estimation for 1-dimensional lists of z-values
#'
#' @description \code{sp.mix.1D} returns LocalFDR estimates and semiparametric mixture density estimates for given 1-dimensional lists of z-values, which are the probit-transformed p-values.
#' For the hypothesis testing \code{sp.mix.1D} uses a two-component semiparametric mixture model to estimate the LocalFDR from the p-values. The two pillars of the proposed approach are Efron’s empirical null principle and log-concave density estimation for the alternative distribution.
#' Estimated prior, null and alternative probability for given points and estimated LocalFDR are returned.
#'
#' @param z Vector which each element indicates z-values, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference of current and previous gamma value is smaller than tol, i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step, then optimization stops. (default: 5e-6)
#' @param max.iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param doplot Boolean parameter that draws histogram and fitted lines for estimated density, if TRUE. (default: TRUE)
#' @param thre.localFDR Threshold of LocalFDR which is used for calculation of false positive rate (FPR) and sensitivity with the threshold set to be \eqn{fdr \leq thre.localFDR}. (default: 0.2)
#'
#' @return Parametrization of f(x) in terms of hyperplanes and function
#'   evaluations y = log(f(x)) \item{aOpt, bOpt}{Analytically normalized
#'   parameters of f(x).} \item{logLike}{Log-likelihood of f(x)} \item{y}{Vector
#'   with values y_i = log(f(X_)) of the normalized density (\eqn{logLike =
#'   \sum(y_i)}).} \item{aOptSparse, bOptSparse}{Sparse parametrization
#'   normalized on the integration grid.}
#' @export
sp.mix.1D <- function(z, tol = 5.0e-6, max.iter = 30, doplot = TRUE, thre.localFDR = 0.2)
{
  #library(LogConcDEAD)
  library(logcondens)

  z <- as.numeric(z)
  n <- length(z)

  ## Initial step
  q0 <- quantile(z, probs = .9)
  p.0 <- mean(z <= q0)
  mu.0 <- mean(z[z <= q0])
  sig.0 <- sd(z[z <= q0])

  mu.1 <- mean(z[z > q0])
  sig.1 <- sd(z[z > q0])
  f1.tilde <- dnorm(z, mean = mu.1, sd = sig.1)

  f <- gam <- rep(0, n)

  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3) | ((k < max.iter) & (!converged)) ) {
    k <- k + 1
    ## E-step
    tmp <- p.0*dnorm(z, mu.0, sig.0)
    new.f <- tmp + (1 - p.0)*f1.tilde
    new.gam <- tmp/new.f

    ## M-step
    w.gam <- new.gam/sum(new.gam, na.rm = TRUE)
    new.mu.0 <- sum(w.gam*z, na.rm = TRUE)
    new.sig.0 <- sqrt(sum(w.gam*(z-new.mu.0)^2, na.rm = TRUE))
    new.p.0 <- mean(new.gam, na.rm = TRUE)

    new.f1.tilde <- rep(0, n)
    which.z <- new.gam <= .95
    weight <- 1 - new.gam[which.z]
    weight <- weight/sum(weight)
    lcd <- activeSetLogCon(x = z[which.z], w = weight)
    new.f1.tilde[which.z] <- exp(lcd$phi)[rank(z[which.z])]
    which.gam <- (new.gam <= 0.9)*(new.gam >= 0.01)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in 1dfdr fit = ", round(diff, 5), "\n")
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde; gam <- new.gam; f <- new.f
  }

  which.z <- gam <= thre.localFDR
  thre <- min(z[which.z])

  if (doplot) {
    hist(z,
         nclass = max(round(length(z)/20), 24),
         probability = TRUE,
         col = "gray", border = "white",
         xlab = "",
         main = "",
         sub = substitute(
           paste(p[0], " = ", p0, ", ",
                 mu[0], " = ", mu0, ", ",
                 sigma[0], " = ", sigma0, ", ",
                 "threshold = ", threshold,
                 sep = ""),
           list(p0 = round(p.0, 2),
                mu0 = round(mu.0, digits = 2),
                sigma0 = round(sig.0, digits = 2),
                threshold = round(thre, digits = 2))))
    rug(z, col = "gray")
    rug(z[which.z], col = 2)
    zs <- sort(z)
    lines(zs, p.0*dnorm(zs, mean = mu.0, sd = sig.0), col = 3, lwd = 2)
    lines(zs, (1-p.0)*f1.tilde[order(z)], col = 2, lwd = 2)
    points(thre, 0, bg = "yellow", col = 2, pch = 25)
  }

  res <- list(p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0, f = f,
              localfdr = gam, iter = k)

  return(res)
}
