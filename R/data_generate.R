library(LaplacesDemon)
library(foreach)
library(doParallel)
library(Matrix)
library(sparseMVN)
library(Rcpp)
sourceCpp("srcs/data_generate.cpp")

spatialMatrix.generate.indices <- function(s, d = 3)
{
  i <- 1:(s ** d - s ** (d - 1) ** 1)
  is <- c(i, i + s ** (d - 1))
  js <- c(i + s ** (d - 1), i)
  m <- rep(1, s); m[1] <- m[s] <- 0
  if(d > 1)
  {
    m <- c()
    small <- spatialMatrix.generate.indices(s, d - 1)
    is.small <- small$is; js.small <- small$js; m.small <- small$m
    for(k in 1:s)
    {
      is <- c(is, is.small + (k - 1) * s ** (d - 1))
      js <- c(js, js.small + (k - 1) * s ** (d - 1))
      if(k == 1 || k == s) m <- c(m, m.small)
      if(k > 1 && k < s) m <- c(m, m.small + 1)
    }
  }
  return(list(is = is, js = js, m = m + 1))
}

spatialMatrix.generate <- function(s, d, sparse = TRUE)
{
  ij <- spatialMatrix.generate.indices(s, d)
  is <- ij$is; js <- ij$js; m <- ij$m
  if(sparse)
  {
    W <- sparseMatrix(is, js, x = 1)
  }
  if(!sparse)
  {
    W <- matrix(0, s ** d)
    for(i in is)
    {
      for(j in js)
      {
        W[i, j] <- 1
      }
    }
  }
  return(list(W = W, M = m))
}

voxels.initialize <- function(s, d, center, act.mean = 3, act.sd = 3, deact.mean = 0, deact.sd = 0.06 ** 0.5)
{
  if(d == 2)
  {
    S <- s ** d
    i <- 0
    act <- rep(FALSE, S); beta <- numeric(S)
    for(s1 in 1:s)
    {
      for(s2 in 1:s)
      {
        i <- i + 1
        if(max(abs(c(s1, s2) - center)) < 4){
          act[i] <- TRUE
          beta[i] <- abs(rnorm(1, 0, act.sd)) + act.mean
        }else{
          beta[i] <- rnorm(1, deact.mean, deact.sd)
        }
      }
    }
  }
  if(d == 3)
  {
    S <- s ** d
    i <- 0
    act <- rep(FALSE, S); beta <- numeric(S)
    for(s1 in 1:s)
    {
      for(s2 in 1:s)
      {
        for(s3 in 1:s)
        {
          i <- i + 1
          if(max(abs(c(s1, s2, s3) - center)) < 4){
            act[i] <- TRUE
            beta[i] <- abs(rnorm(1, 0, act.sd)) + act.mean
          }else{
            beta[i] <- rnorm(1, deact.mean, deact.sd)
          }
        }
      }
    }
  }
  return(list(beta = beta, act = act))
}

make.M <- function(s, d, begind)
{
  if(d == 0) return(begind)
  a <- make.M(s, d - 1, begind + 1)
  M <- rep(a, s)
  M[1:s ** (d - 1)] <- M[1:s ** (d - 1)] - 1
  M[1:s ** (d - 1) + s ** (d - 1) * (s - 1)] <- M[1:s ** (d - 1) + s ** (d - 1) * (s - 1)] - 1
  return(M)
}
find.mu <- function(s, d, delta, begind, x)
{
  if(d == 0) return(x/begind)
  t1 <- Sys.time()
  W <- as.matrix(spatialMatrix.generate(s, d)$W)
  M <- make.M(s, d, begind)
  Sigma <- solve(diag(M) - delta * W)
  tt <<- tt + (Sys.time() - t1); t1 <- Sys.time()
  y <- x %*% Sigma
  tt1 <<- tt1 + (Sys.time() - t1)
  return(y)
  # solve(a, b, sparse = FALSE, tol = .Machine$double.eps, ...)
}
sampling <- function(n, s, d = 2, delta = 0.3, begind = d)
{
  if(d == 0) return(matrix(rnorm(n), n, 1) / sqrt(begind))
  res <- matrix(0, n, s ** d)
  ind <- 1:s ** (d - 1)

  boundary <- sampling(2 * n, s, d - 1, delta, begind)
  nonboundary <- sampling((s - 2) * n, s, d - 1, delta, begind + 1)

  for(i in 1:n)
  {
    res[i, ind] <- boundary[i, ]
    for(j in 2:(s - 1)) res[i, (j - 1) * s ** (d - 1) + ind] <- nonboundary[(i - 1) * (s - 2) + j - 1, ]
    res[i, ind + (s - 1) * s ** (d - 1)] <- boundary[i + n, ]
  }

  if(d == 1)
  {
    for(j in 2:s)
    {
      i1 <- (j - 2) * s ** (d - 1) + ind
      i2 <- (j - 1) * s ** (d - 1) + ind
      res[, i2] <- res[, i2] +  delta * find.mu(s, d - 1, delta, begind + (j != s), res[, i1])
    }
    return(res)
  }

  W <- as.matrix(spatialMatrix.generate(s, d - 1)$W)
  M <- make.M(s, d - 1, begind + 1)
  Omega.now <- diag(M) - delta * W
  # solve(a, b, sparse = FALSE, tol = .Machine$double.eps, ...)

  for(j in 2:s)
  {
    if(j == s) Omega.now <- Omega.now - diag(s ** (d - 1))
    for(i in 1:n)
    {
      i1 <- (j - 2) * s ** (d - 1) + ind
      i2 <- (j - 1) * s ** (d - 1) + ind
      # res[i, i2] <- res[i, i2] +  delta * find.mu(s, d - 1, delta, begind + (j != s), res[i, i1])
      res[i, i2] <- res[i, i2] + solve(Omega.now, res[i, i1], sparse = TRUE)
    }
  }

  return(res)
}

data_generate <- function(theta, covariates, N, s, d, T, setting = 0)
{
  epsilon <- t(sampling(T, s, d, theta$delta))
  return(array(data_generate_cpp(theta, covariates, N, s ** d, T, 3, epsilon), dim = c(N, S, T)))
}
