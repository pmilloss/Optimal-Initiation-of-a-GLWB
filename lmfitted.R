lm.fitted <- function (x, y) 
{
  be <- solve(crossprod(x), crossprod(x, y))
  ff <- x %*% be
  return(ff)
}