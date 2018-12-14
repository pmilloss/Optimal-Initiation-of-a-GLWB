trapezi <- function(f, d)
  {
  0.5 * d * sum(f[- length(f)] + f[- 1])
  }
  
trapezi.cumulati <- function(f, d)
  {
  cf <- cumsum(f)
  return(0.5 * d * (cf[- 1] + cf[- length(cf)] - f[1]))
  }

confronto <- function(x, y)
  {
  l <- length(x)
  if (x[l] > y) return(min(which(x > y))) else return(l + 1)
  # l+1=survival
  }

weibull.intensity <- function(t, x, theta, c)
  {
  (c / (theta ^ c)) * ((x + t) ^ (c - 1))
  }