XM1.0 <- function(X_old, M_old, dS_S, d, psi, phi, beta)
{
  X_new <- pmax(X_old * (1 + dS_S - psi * d) - phi * d * M_old, 0)
  M_new <- pmax(M_old * (1 + beta * d), X_new)
  
  return(cbind(X_new, M_new))
}

XMt1.0 <- function(Xt_old, Mt_old, dS_S, d, psi, phi, beta = 0, gt)
{
  Xt_new <- pmax(Xt_old * (1 + dS_S - psi * d) - (phi + gt) * d * Mt_old, 0)
  Mt_new <- pmax(Mt_old * (1 + beta * d), Xt_new)
  
  return(cbind(Xt_new, Mt_new))
}