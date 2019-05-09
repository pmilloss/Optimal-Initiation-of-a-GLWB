simula_dS_S_GBM <- function(nsim, N, r, sigma)
{
  dW <- matrix(rnorm(nsim * N), nrow = N, ncol = nsim)
  dS_S <- exp((r - 0.5 * sigma ^ 2) + sigma * dW) - 1
  return(dS_S)  
}