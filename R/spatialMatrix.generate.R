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
