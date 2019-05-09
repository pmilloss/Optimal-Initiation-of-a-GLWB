optimal_initiation_surrender <- function(nsim, X, M, state_var = NULL, G, basis_fun, Nbar, d, tau, B)
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
  
  lambda <- tau
  payoff <- X[cbind(tau + 1, 1 : nsim)]
  # payoff in case of death or initiation
  
  # Nbar <- max(tau)
  
  for (t in (Nbar - 1) : 0)
  {
    ind_surv_t <- which(tau > t)
    # (indices of) survivors
    
    init_value_surr <- G[t + 1, ind_surv_t]
    
    cont_value <- payoff[ind_surv_t] * B[t + 1, ind_surv_t] / B[cbind(lambda[ind_surv_t] + 1, ind_surv_t)]
    
    if(is.null(state_var)) Z <- cbind(1, X[t + 1, ind_surv_t], M[t + 1, ind_surv_t]) else Z <- cbind(1, X[t + 1, ind_surv_t], M[t + 1, ind_surv_t], sapply(state_var, extract_row_col, row = t + 1, col = ind_surv_t))
    
    reg_mat <- t(apply(Z, 1, basis_fun))
    # apply "basis_functions" to each row of Z to build the regression matrix
    
    # LSMC to calculate the continuation value
    cont_value_LSMC <- lm.fit(reg_mat, cont_value)$fitted.values
    
    # LSMC to calculate the initiation value
    init_value_surr_LSMC <- lm.fit(reg_mat, init_value_surr)$fitted.values
    
    index_init <- (init_value_surr_LSMC > cont_value_LSMC)
    
    payoff[ind_surv_t][index_init] <- init_value_surr[index_init]
    lambda[ind_surv_t][index_init] <- t
    
  }
  contract_value_0 <- mean(payoff / B[lambda + 1, ])
  SE <- sd(payoff / B[lambda + 1, ]) / sqrt(nsim)
  return(list(contract_value_0, SE, lambda))
}