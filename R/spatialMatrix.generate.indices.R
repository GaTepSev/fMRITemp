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
