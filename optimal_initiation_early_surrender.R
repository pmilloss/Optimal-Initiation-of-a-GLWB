optimal_initiation_early_surrender <- function(nsim, X, M, state_variables, F, basis_functions, n_basis_functions, N, d, tau, B, penalty)
{
  # calculate contract value at time 0 that can be initiated optimally and early surrendered
  
  # state_variables: should include all state variables but the personal and base account
  lambda <- pai <- tau
  payoff_initiation_surrender <- X[cbind(tau + 1, 1 : nsim)]
  
  Nbar <- max(tau)
  
  for (n in (Nbar - 1) : 0)
  {
    ind_surv <- which(tau > n)

    initiation_value <- F[n + 1, ind_surv]

    continuation_value <- payoff_initiation_surrender[ind_surv] * B[n + 1, ind_surv] / B[cbind(pmin(lambda[ind_surv_j], pai[ind_surv_j]) + 1, ind_surv)]
    
    # regression_matrix[j, ] <- basis_functions(c(X[n + 1, ind_surv_j], M[n + 1, ind_surv_j])) # AGGIUNGERE ALTRE VARIABILI DI 
        
    state_variables<- cbind(X[n + 1, ind_surv], M[n + 1, ind_surv])
    
    regression_matrix <- t(apply(X = state_variables, MARGIN = 1, FUN = "basis_functions"))
    # basis_functions(c(X[n + 1, ind_surv_j], M[n + 1, ind_surv_j])) # AGGIUNGERE ALTRE VARIABILI DI STATO
    
    # initiation_value <- continuation_value <- surrender_value <- c()
    # regression_matrix <- matrix(nrow = l_surv, ncol = n_basis_functions)
    
    # initiation_value <- F[n + 1, ind_surv]
    
    # for(j in 1 : l_surv)
    # {
    #   ind_surv_j <- ind_surv[j]
    #   
    #   continuation_value[j] <- payoff_initiation_surrender[ind_surv_j] * B[n + 1, ind_surv_j] / B[pmin(lambda[ind_surv_j], pai[ind_surv_j]) + 1, ind_surv_j]
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
    
    surrender_value <- X[n + 1, ind_surv] * (1 - penalty)
    
    index_initiation <- (initiation_value_LSMC > pmax(continuation_value_LSMC, surrender_value))

    index_surrender <- (surrender_value > pmax(continuation_value_LSMC, initiation_value_LSMC))    
    # indiciRiscatti <- (valoreRiscatto > valoreContinuazioneStimato)
    payoff_initiation_surrender[ind_surv][index_initiation] <- initiation_value[index_initiation]
    lambda[ind_surv][index_initiation] <- n
  
    payoff_initiation_surrender[ind_surv][index_surrender] <- surrender_value[index_surrender]
    pai[ind_surv][index_surrender] <- n
    }
  
  contract_value_0 <- mean(payoff_initiation_surrender / B[pmin(lambda, pai) + 1, ])
  
  SE <- sd(payoff_initiation_surrender / B[pmin(lambda, pai) + 1, ]) / sqrt(nsim)
    
  return(list(contract_value_0, SE, lambda, pai))
}


system.time(res3 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 1))

system.time(res3 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 0.2))

XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phi = 0.03, beta)

FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phi = 0.03, beta = 0, gt = g, tau, B)


system.time(res4 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 0))

system.time(res5 <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 0.04))


# 
# system.time(res2 <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B))


