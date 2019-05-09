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

optimal_initiation_early_surrender <- function(nsim, X, M, state_var = NULL, FF, basis_fun, tau, Nbar, B, penalty)
  {
  penalty <- rep(penalty, length.out = Nbar)
  lambda <- pai <- tau
  payoff <- X[cbind(tau + 1, 1 : nsim)]
  
  for (t in (Nbar - 1) : 0)
  {
    ind_surv_t <- which(tau > t)
    init_value <- FF[t + 1, ind_surv_t]
    
    cont_value <- payoff[ind_surv_t] * B[t + 1, ind_surv_t] / B[cbind(pmin(lambda[ind_surv_t], pai[ind_surv_t]) + 1, ind_surv_t)]
    
    if(is.null(state_var)) Z <- cbind(1, X[t + 1, ind_surv_t], M[t + 1, ind_surv_t]) else Z <- cbind(1, X[t + 1, ind_surv_t], M[t + 1, ind_surv_t], sapply(state_var, extract_row_col, row = t + 1, col = ind_surv_t))

    reg_mat <- t(apply(Z, 1, basis_fun))
    cont_value_LSMC <- lm.fit(reg_mat, cont_value)$fitted.values
    init_value_LSMC <- lm.fit(reg_mat, init_value)$fitted.values
    surr_value <- X[t + 1, ind_surv_t] * (1 - penalty[t + 1])
    
    index_init <- (init_value_LSMC > pmax(cont_value_LSMC, surr_value))
    index_surr <- (surr_value > pmax(cont_value_LSMC, init_value_LSMC))
    
    payoff[ind_surv_t][index_init] <- init_value[index_init]
    lambda[ind_surv_t][index_init] <- t
    payoff[ind_surv_t][index_surr] <- surr_value[index_surr]
    pai[ind_surv_t][index_surr] <- t
    }
  contract_value_0 <- mean(payoff / B[pmin(lambda, pai) + 1, ])
  SE <- sd(payoff / B[pmin(lambda, pai) + 1, ]) / sqrt(nsim)
  return(list(contract_value_0, SE, lambda, pai))
}