# calculate simulated values of X and M given simulated returns dS_S

# nsim: number of simulations, a scalar
# X0: initial value of the personal account, a scalar
# M0: initial value of the base account, a scalar
# dS_S: a (N x nsim) - matrix with the simulated log returns of the reference fund 
# N: the number of periods, a scalar (integer)
# d: the time step (a fraction of year)
# psi: the management fee
# phi: the insurance fee
# beta: the roll-up rate

simulate_XM <- function(nsim, X0 = 100, M0 = X0, dS_S, Nbar, psi, phi, beta)
{
  X <- M <- matrix(0, nrow = Nbar + 1, ncol = nsim)
  rownames(X) <- rownames(M) <- seq(0, Nbar)
  # set up matrices for X and M
  
  X[1, ] <- X0
  M[1, ] <- M0
  # initialize first row (time 0)
  
  for(t in 1 : Nbar)
  {
    # given X[t, ], M[t, ] (time t - 1) calculate X[t + 1, ], M[t + 1] (time t)
    
    X[t + 1, ] <- pmax(X[t, ] * (1 + dS_S[t, ] - psi) - phi * M[t, ], 0)
    M[t + 1, ] <- pmax(M[t, ] * (1 + beta), X[t + 1, ])
  }
  return(list(X, M))
}
