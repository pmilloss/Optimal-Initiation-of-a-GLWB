simulate_mu_tau <- function(nsim, N, M, k_mu, sigma_mu, x0, b, m)
{
  l <- N * M
  delta <- 1 / M
  
  mt <- function(t)exp((x0 + t - m) / b) / b
  mu0 <- mt(0)
  mu <- matrix(nrow = N, ncol = nsim)
  
  tau <- mu_int <- vector(length = nsim, mode = 'numeric')
  exp_unit <- rexp(nsim)
  mu_temp <- matrix(nrow = M + 1, ncol = nsim)
  mu[1, ] <- mu_temp[1, ] <- mu0
  mu_int <- 0
  
  for (t in 1 : N)
  {
    for (u in 1 : M)
    {
      dW <- rnorm(nsim, sd = sqrt(delta))
      mu_temp[u + 1, ] <- CIR1(mu_temp[u, ], k_mu, mt((t - 1) + (u - 1) * delta), sigma_mu, delta, dW)
    }
    mu[t, ] <- mu_temp[M + 1, ]
    mu_int_temp <- apply(mu_temp, 2, trap, d = delta)
    tau[(mu_int <= exp_unit) & (mu_int + mu_int_temp > exp_unit)] <- t
    mu_temp[1, ] <- mu[t, ]
    mu_int <- mu_int + mu_int_temp 
  }
  tau[mu_int <= exp_unit] <- N
  return(list(mu, tau))
}


# simulazione <- function(parm_sim, parm_contratto, parm_mu, parm_r, parm_V, parm_S, parm_A, odecoef)
# {
#   k_mu <- parm_mu[1]; sigma_mu <- parm_mu[2]; x <- parm_mu[3];
#   c <- parm_mu[4]; theta <- parm_mu[5]
#   r0 <- parm_r[1]; k_r <- parm_r[2]; theta_r <- parm_r[3]; sigma_r <- parm_r[4]
#   V0 <- parm_V[1]; k_V <- parm_V[2]; theta_V <- parm_V[3]; sigma_V <- parm_V[4]
#   S0 <- parm_S[1]; rho_VS <- parm_S[2]
#   A0 <- parm_A
#   
#   l <- N*M
#   Delta <- maturita/N; delta <- Delta/M
#   
#   m <- function(t)weibull.intensity(t, x = x, theta = theta, c = c); mu0 <- m(0)
#   kk <- 1-phi*delta
#   
#   xx <- seq(0, by = delta, length.out = M+1)
#   ind <- seq(M, l, by = M)
#   indxx <- c(1, ind[-N]+1)
#   
#   mu <- r <- fattoriSconto <- V <- A <- integrale_mu <- integrale_r <- matrix(nrow = N, ncol = numSim)
#   
#   payoff <- decessi <- vector(length = numSim, mode = 'numeric')
#   expUnit <- rexp(numSim)
#   mu_temp <- r_temp <- V_temp <- A_temp <- matrix(nrow = l+1, ncol = numSim)
#   mu_temp[1, ] <- mu0; r_temp[1, ] <- r0; V_temp[1, ] <- V0; A_temp[1, ] <- A0
#   
#   for (j in 1:l)
#   {
#     dW <- matrix(rnorm(4*(numSim/2), sd = sqrt(delta)), ncol = numSim/2, nrow = 4)
#     dW <- cbind(dW, -dW)
#     mu_temp[j+1, ] <- CIR1(mu_temp[j, ], k_mu, m((j-1)*delta), sigma_mu, delta, dW[1, ])
#     r_temp[j+1, ] <- CIR1(r_temp[j, ], k_r, theta_r, sigma_r, delta, dW[2, ])
#     V_temp[j+1, ] <- CIR1(V_temp[j, ], k_V, theta_V, sigma_V, delta, dW[3, ])
#     dS.S <- dS1.S0(r_temp[j, ], V_temp[j, ], delta, dW[3, ], dW[4, ], rho_VS)
#     A_temp[j+1, ] <- A_temp[j, ]*(kk+dS.S)
#   }
#   mu <- mu_temp[ind+1, ]; r <- r_temp[ind+1, ]; V <- V_temp[ind+1, ]; A <- A_temp[ind+1, ]
#   for (h in 1:N)
#   {
#     integrale_mu[h, ] <- apply(mu_temp[seq(indxx[h], length.out = M+1), ], 2, trapezi, d = delta)
#     integrale_r[h, ] <- apply(r_temp[seq(indxx[h], length.out = M+1), ], 2, trapezi, d = delta)
#   }
#   integrale_mu <- apply(integrale_mu, 2, cumsum)
#   for (n in 1:numSim)decessi[n] <- confronto(integrale_mu[, n], expUnit[n])
#   fattoriSconto <- exp(-integrale_r)
#   
#   ind.sopravvive <- which(decessi > N); decessi[ind.sopravvive] <- N
#   payoff <- pmax(A[cbind(decessi, 1:numSim)], minGarantito[decessi])
#   payoff[ind.sopravvive] <- A[N, ind.sopravvive]
#   return(list(mu, decessi, r, fattoriSconto, V, A, payoff))
# }