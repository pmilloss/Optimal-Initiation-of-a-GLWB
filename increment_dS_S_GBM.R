simula_dS_S_GBM <- function(nsim, N, r, sigma, d)
{
  dW <- matrix(rnorm(nsim * N, 0, sqrt(d)), nrow = N, ncol = nsim)
  dS_S <- exp((r - 0.5 * sigma ^ 2) * d + sigma * dW) - 1
  return(dS_S)  
}