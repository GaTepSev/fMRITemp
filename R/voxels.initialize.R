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
