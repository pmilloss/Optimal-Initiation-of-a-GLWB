optimal_initiation_early_surrender <- function(nsim, X, M, state_var = NULL, F, basis_fun, n_basis_fun, N, d, tau, B, penalty)
{
  # calculate contract value at time 0 that can be initiated optimally and early surrendered
  
  # state_variables: should include all state variables but the personal and base account
  
  # nsim: number of simulations
  # X: a (N+1) x nsim matrix with personal account simulations
  # M: a (N+1) x nsim matrix with base amount simulations,
  # state_variables: a list of (N+1) x nsim matrices with additional state variables, or NULL
  # F: a N x nsim matrix with present values of cashflows received in case of initiation
  # basis_functions: a vector function returning a vector of basis functions
  # n_basis_functions: the number of basis functions
  # N: the number of periods spanning the remaining lifetime of the policyholder
  # d: the interval length
  # tau: a vector of integer based death times
  # B: an (N+1) x nsim matrix with money market instrument values
  # penalty: a vector with penalty

  penalty <- rep(penalty, length.out = N * d)
  lambda <- pai <- tau
  payoff <- X[cbind(tau + 1, 1 : nsim)]
  # payoff in case of death, initiation or surrender
  
  Nbar <- max(tau)
  
  for (n in (Nbar - 1) : 0)
  {
    ind_surv <- which(tau > n)
    # (indices of) survivors

    init_value <- F[n + 1, ind_surv]

    cont_value <- payoff[ind_surv] * B[n + 1, ind_surv] / B[cbind(pmin(lambda[ind_surv], pai[ind_surv]) + 1, ind_surv)]
    
    # regression_matrix[j, ] <- basis_functions(c(X[n + 1, ind_surv_j], M[n + 1, ind_surv_j])) # AGGIUNGERE ALTRE VARIABILI DI STATO
    Z <- cbind(X[n + 1, ind_surv], M[n + 1, ind_surv])
    # matrix of state variables
    
    reg_mat <- t(apply(Z, 1, basis_fun))
    # apply "basis_functions" to each row of Z to build the regression matrix
    
    # LSMC to calculate the continuation value
    cont_value_LSMC <- lm.fit(reg_mat, cont_value)$fitted.values
    
    # LSMC to calculate the initiation value
    init_value_LSMC <- lm.fit(reg_mat, init_value)$fitted.values
    
    surr_value <- X[n + 1, ind_surv] * (1 - penalty[n + 1])
    
    index_init <- (init_value_LSMC > pmax(cont_value_LSMC, surr_value))

    index_surr <- (surr_value > pmax(cont_value_LSMC, init_value_LSMC))    
    payoff[ind_surv][index_init] <- init_value[index_init]
    lambda[ind_surv][index_init] <- n
  
    payoff[ind_surv][index_surr] <- surr_value[index_surr]
    pai[ind_surv][index_surr] <- n
    }
  
  contract_value_0 <- mean(payoff / B[pmin(lambda, pai) + 1, ])
  
  SE <- sd(payoff / B[pmin(lambda, pai) + 1, ]) / sqrt(nsim)
    
  return(list(contract_value_0, SE, lambda, pai))
}


# system.time(res3 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 1))
# 
# system.time(res3 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 0.2))
# 
# XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phi = 0.03, beta)
# 
# FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phi = 0.03, beta = 0, gt = g, tau, B)
# 
# 
# system.time(res4 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 0))
# 
# system.time(res5 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 0.04))


# 
# system.time(res2 <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B))


