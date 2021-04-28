# random variable generator for Dirichlet negative multinomial distribution
require(numDeriv)
rDirNegMulti <- function(n, nSuccesses = 1, dirPars = c(1, 1)) {
  m <- length(dirPars) - 1
  lambda0 <- rgamma(n, shape = dirPars[1], rate = nSuccesses)
  lambdaj <- matrix(
    t(sapply(lambda0, function(rate) rgamma(m, shape = dirPars[-1], rate = rate))),
    nrow = n,
    ncol = m
  )
  t <- rgamma(n, shape = nSuccesses, rate = nSuccesses)
  y <- matrix(
    rpois(n * m, lambda = lambdaj * t),
    nrow = n,
    ncol = m
    )
    y 
}

#log-likelihood function
dDirNegMulti <- function(x, nSuccesses = 1, dirPars = c(1, 1), log = FALSE) {

  x <- cbind(nSuccesses, x)
  dirParsMat <- matrix(dirPars, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  xDirParsMat <- x + dirParsMat
  ld <- lgamma(apply(x, 1, sum)) - lgamma(x[, 1]) -
    apply(lgamma(x[, -1, drop = FALSE] + 1), 1, sum)
  ld <- ld + lgamma(sum(dirPars)) - sum(lgamma(dirPars))
  ld <- ld + apply(lgamma(xDirParsMat), 1, sum) -
    lgamma(apply(xDirParsMat, 1, sum))
  if(log) ld else exp(ld)
}

# Jacobian function
jDirNegMulti <- function(x, nSuccesses = 1, dirPars = rep(1, m + 1)) {

  n <- nrow(x)
  m <- ncol(x)
  x <- cbind(nSuccesses, x)
  dirParsMat <- matrix(dirPars, nrow = n, ncol = m + 1, byrow = TRUE)
  xDirParsMat <- x + dirParsMat
  cbind(
    nSuccesses = digamma(apply(x, 1, sum)) - digamma(x[, 1]) +
      digamma(xDirParsMat[, 1]) - digamma(apply(xDirParsMat, 1, sum)),
    dirPars = digamma(apply(dirParsMat, 1, sum)) -
      digamma(apply(xDirParsMat, 1, sum)) +
      matrix(sapply(1:(m + 1), function(j) {
        digamma(xDirParsMat[, j]) - digamma(dirParsMat[, j])
      }), ncol = m + 1)
  )
}

## fit, using regression parametrisation
glmDirNegMulti <- function(formula, group, data = NULL, family = poisson(), ...) {
  modFormula <- formula
  modFormula[[3]] <- as.call(list(as.name("+"), modFormula[[3]], as.name(group)))
  mf <- model.frame(modFormula, data)
  glmFit <- glm(formula = formula, data = mf, family = family)
  logLik <- function(par, unconstrPars = FALSE) {
    if(unconstrPars) {
      phiNM <- exp(par[1])
      phiNB <- exp(par[2])
    } else {
      phiNM <- par[1]
      phiNB <- par[2]
    }
    ll <- by(mf, mf[[group]], function(sf) {
      mu <- family$linkinv(model.matrix(formula, sf) %*% par[-(1:2)])
      y <- model.response(sf)
      dDirNegMulti(matrix(y, nrow = 1), nSuccesses = phiNM^(-1),
                   dirPars = c(phiNM^(-1) * phiNB^(-1) + 1,
                               mu * phiNB^(-1)),
                   log = TRUE)
    })
    sum(ll)
  }
  jacob <- function(par, unconstrPars = FALSE) {
    if(unconstrPars) {
      phiNM <- exp(par[1])
      phiNB <- exp(par[2])
    } else {
      phiNM <- par[1]
      phiNB <- par[2]
    }
    p <- length(coef(glmFit))
    ja <- by(mf, mf[[group]], function(sf) {
      mu <- family$linkinv(model.matrix(formula, sf) %*% par[-(1:2)])
      y <- model.response(sf)
      m <- length(y)
      regPars.unconstrPars <- diag(c(phiNM, phiNB, rep(1, p)))
      meanPars.regPars <- matrix(0, nrow = p + 2, ncol = m + 2)
      meanPars.regPars[cbind(1:2, 1:2)] <- -c(phiNM, phiNB)^(-2)
      meanPars.regPars[-(1:2), -(1:2)] <- t(as.vector(
        family$mu.eta(family$linkfun(mu))) *
          model.matrix(formula, sf))
      natPars.meanPars <- diag(c(1, phiNM^(-1), rep(phiNB^(-1), m)))
      natPars.meanPars[1, 2] <- phiNB^(-1)
      natPars.meanPars[2, -(1:2)] <- mu
      l.natPars <- t(jDirNegMulti(matrix(y, nrow = 1),
                                  nSuccesses = phiNM^(-1),
                                  dirPars = c(phiNM^(-1) * phiNB^(-1) + 1,
                                              mu *phiNB^(-1))))
      if(unconstrPars) {
        regPars.unconstrPars %*% meanPars.regPars %*%
          natPars.meanPars %*% l.natPars
      } else {
        meanPars.regPars %*% natPars.meanPars %*% l.natPars
      }
    })
    apply(do.call(cbind, ja), 1, sum)
  }
  constrFit <- optim(
    par = c(0.5, 0.5, coef(glmFit)),
    fn = logLik,
    gr = jacob,
    unconstrPars = FALSE,
    method = "L-BFGS-B",
    lower = c(1e-2, 1e-2, rep(-Inf, length(coef(glmFit)))),
    control = list(fnscale = -1)
  )
  unconstrMLE <- c(log(constrFit$par[1:2]), constrFit$par[-(1:2)])
  hess <- hessian(logLik, unconstrMLE, unconstrPars = TRUE)
  vc <- chol2inv(chol(-hess))
  glmVC <- vc[-(1:2), -(1:2)]
  rownames(glmVC) <- colnames(glmVC) <- names(coef(glmFit))
  phiNM <- exp(unconstrMLE[1])
  attr(phiNM, "CI95") <- exp(unconstrMLE[1] +
                               c(-1, 1) * qnorm(0.975) * sqrt(vc[1, 1]))
  phiNB <- exp(unconstrMLE[2])
  attr(phiNB, "CI95") <- exp(unconstrMLE[2] +
                               c(-1, 1) * qnorm(0.975) * sqrt(vc[2, 2]))
  list(phiNM = phiNM,
       phiNB = phiNB,
       coef = unconstrMLE[-(1:2)],
       vcov = glmVC,
       logLik = constrFit$value,
       convergence = constrFit$message)
}



#an example
set.seed(1)
fr <- data.frame(i = rep(1:22, times = c(rep(3:4, times = c(12, 10)))))
fr <- within(fr, {
  x1 <- rbinom(22, size = 1, prob = 0.5)[i]
  x2 <- rbinom(22, size = 1, prob = 0.5)[i]
  x3 <- rbinom(length(i), size = 1, prob = 0.5)
})
fr <- do.call(rbind, by(fr, fr$i, function(sf) {
  m <- nrow(sf)
  mu <- exp(1 + 0.93 * sf$x1 + 0.68 * sf$x2 + 0.36 * sf$x3)
  y <- rDirNegMulti(n = 1, nSuccesses = 16.09, dirPars = c(16.09 * 0.357 + 1, mu * 0.357))
  sf$y <- as.numeric(y)
  sf
}))
fit <- glmDirNegMulti(y ~ x1 + x2 + x3, group = "i", data = fr)
print(fit)