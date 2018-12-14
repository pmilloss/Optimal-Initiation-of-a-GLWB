# optimal_initiation <- function(mu, decessi, r, fattoriSconto, V, A, payoff, funzioniBase, numFunzioniBase, parm_sim, penalita)
optimal_initiation <- function(nsim, X, M, state_variables, F, basis_functions, n_basis_functions, N, d, tau, B)
{
  # calculate contract value at time 0 that can be initiated optimally but not surrendered
  
  # state_variables: should include all state variables but the personal and base account
  lambda <- tau
  payoff_initiation <- X[cbind(tau + 1, 1 : nsim)]
  
  Nbar <- max(tau)
  
  for (n in (Nbar - 1) : 0)
  {
    ind_surv <- which(tau > n)
    # l_surv <- length(ind_surv)
    # initiation_value <- continuation_value <- c()
    # regression_matrix <- matrix(nrow = l_surv, ncol = n_basis_functions)
    
    initiation_value <- F[n + 1, ind_surv]

    continuation_value <- payoff_initiation[ind_surv] * B[n + 1, ind_surv] / B[cbind(lambda[ind_surv] + 1, ind_surv)]

    state_variables<- cbind(X[n + 1, ind_surv], M[n + 1, ind_surv])
    
    regression_matrix <- t(apply(X = state_variables, MARGIN = 1, FUN = "basis_functions"))
      # basis_functions(c(X[n + 1, ind_surv_j], M[n + 1, ind_surv_j])) # AGGIUNGERE ALTRE VARIABILI DI STATO
    
        
    # for(j in 1 : l_surv)
    # {
    #   ind_surv_j <- ind_surv[j]
    #   
    #   continuation_value[j] <- payoff_initiation[ind_surv_j] * B[n + 1, ind_surv_j] / B[lambda[ind_surv_j] + 1, ind_surv_j]
    #   
    #   regression_matrix[j, ] <- basis_functions(c(X[n + 1, ind_surv_j], M[n + 1, ind_surv_j])) # AGGIUNGERE ALTRE VARIABILI DI STATO
    # }
    
    # continuation_value_LSMC <- as.numeric(fitted(lm(continuation_value ~ regression_matrix)))
    # 
    # initiation_value_LSMC <- as.numeric(fitted(lm(initiation_value ~ regression_matrix)))
    
    continuation_value_LSMC <- lm.fit(regression_matrix, continuation_value)$fitted.values

    initiation_value_LSMC <- lm.fit(regression_matrix, initiation_value)$fitted.values
    
    # continuation_value_LSMC <- lm.fitted(regression_matrix, continuation_value)
    # 
    # initiation_value_LSMC <- lm.fitted(regression_matrix, initiation_value)
    
    index_initiation <- (initiation_value_LSMC > continuation_value_LSMC)
    
    # indiciRiscatti <- (valoreRiscatto > valoreContinuazioneStimato)
    payoff_initiation[ind_surv][index_initiation] <- initiation_value[index_initiation]
    lambda[ind_surv][index_initiation] <- n
  }
  
  contract_value_0 <- mean(payoff_initiation / B[lambda + 1, ])
  
  SE <- sd(payoff_initiation / B[lambda + 1, ]) / sqrt(nsim)
  
  return(list(contract_value_0, SE, lambda))
}


res1 <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B)
   
system.time(res2 <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B))


