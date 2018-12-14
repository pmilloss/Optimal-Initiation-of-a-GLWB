# calculate fair insurance fees in the initiation only and initiation + early surrender case

X0 <- 1
M0 <- X0

r <- 0.03
sigma <- 0.20

psi <- 0.0075
phi <- 0

g <- 0.05

beta <- 0.05

nsim <- 10000
N <- 60
d <- 1


set.seed(0)

B <- matrix(exp(r * (0 : N)), nrow = N + 1, ncol = nsim, byrow = FALSE)

dS_S <- simula_dS_S_GBM(nsim, N, r, sigma, d)

x0 <- 60; m <- 87.25; b <- 9.5
tau <- simulate_tau_gompertz(nsim = nsim, x0 = x0, b = b, m = m)
tau <- ceiling(tau)





fff <- function(){

  fairphi <- c()
  

  phiV <- c(0, 0.1, seq(0.01, 0.09, by = 0.01))
  VV <- vector(mode = "numeric", length = length(phiV))
    
    for (j in 1 : length(phiV))
    {
      
      XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phiV[j], beta)
      
      FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phiV[j], beta = 0, gt = g, tau, B)
      
      VV[j] <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B)[[1]]
      
      print(c(j, VV[j]))
      
      if (j == 1 & VV[1] < X0) {fairphi <- -99; break()}
      if (j == 2 & VV[2] > X0) {fairphi <- 99; break()}
      
      
    }
    if (fairphi %in% c(-99, 99)) next()
    
    splineV <- function(fff)splinefun(x = phiV, y = VV, method = "hyman" )(fff) - X0
    fairphi <- uniroot(f = splineV, interval = range(phiV))$root
    print(fairphi)
  return(fairphi)
}










fff <- function(){
betaV <- seq(0, 0.1, by = 0.02)
# gV <- seq(0.01, 0.1, by = 0.01)  
  
fairphi <- vector(mode = "numeric", length = length(betaV))
  
for (i in 1 : length(betaV))
{
  
  phiV <- c(0, 0.1, seq(0.01, 0.09, by = 0.01))
  VV <- vector(mode = "numeric", length = length(phiV))
    
  for (j in 1 : length(phiV))
  {

    XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phiV[j], betaV[i])
  
    FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phiV[j], beta = 0, gt = g, tau, B)

    VV[j] <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B)[[1]]

    print(c(i, j, VV[j]))
    
    if (j == 1 & VV[1] < X0) {fairphi[i] <- -99; break()}
    if (j == 2 & VV[2] > X0) {fairphi[i] <- 99; break()}

    
  }
  if (fairphi[i] %in% c(-99, 99)) next()
  
  splineV <- function(fff)splinefun(x = phiV, y = VV)(fff) - X0
  fairphi[i] <- uniroot(f = splineV, interval = range(phiV))$root
  print(fairphi[i])
}
  return(fairphi)
}




fff2 <- function(){
  gV <- seq(0.04, 0.06, by = 0.005)  
  
  fairphi <- vector(mode = "numeric", length = length(gV))
  
  for (i in 1 : length(gV))
  {
    
    phiV <- c(0, 0.5, seq(0.001, 0.01, by = 0.001), 0.05, 0.1)
    VV <- vector(mode = "numeric", length = length(phiV))
    
    for (j in 1 : length(phiV))
    {
      
      XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phiV[j], beta)
      
      FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phiV[j], beta = 0, gt = gV[i], tau, B)
      
      VV[j] <- optimal_initiation(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B)[[1]]
      
      print(c(gV[i], phiV[j], VV[j]))
      
      if (j == 1 & VV[1] < X0) {fairphi[i] <- -99; break()}
      if (j == 2 & VV[2] > X0) {fairphi[i] <- 99; break()}
      
    }
    if (fairphi[i] %in% c(-99, 99)) next()
    
    # splineV <- function(fff)splinefun(x = phiV, y = VV)(fff) - X0
    splineV <- function(fff)approxfun(x = phiV, y = VV)(fff) - X0
    fairphi[i] <- uniroot(f = splineV, interval = range(phiV))$root
    print(fairphi[i])
  }
  return(list(gV, fairphi))
}



fff2s <- function(){
  gV <- seq(0.04, 0.06, by = 0.005)  
  
  fairphi <- vector(mode = "numeric", length = length(gV))
  
  for (i in 1 : length(gV))
  {
    
    phiV <- c(0, 0.5, seq(0.001, 0.01, by = 0.001), 0.05, 0.1)
    VV <- vector(mode = "numeric", length = length(phiV))
    
    for (j in 1 : length(phiV))
    {
      
      XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phiV[j], beta)
      
      FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phiV[j], beta = 0, gt = gV[i], tau, B)
      
      VV[j] <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 1)[[1]]
      
      print(c(gV[i], phiV[j], VV[j]))
      
      if (j == 1 & VV[1] < X0) {fairphi[i] <- -99; break()}
      if (j == 2 & VV[2] > X0) {fairphi[i] <- 99; break()}
      
    }
    if (fairphi[i] %in% c(-99, 99)) next()
    
    # splineV <- function(fff)splinefun(x = phiV, y = VV)(fff) - X0
    splineV <- function(fff)approxfun(x = phiV, y = VV)(fff) - X0
    fairphi[i] <- uniroot(f = splineV, interval = range(phiV))$root
    print(fairphi[i])
  }
  return(list(gV, fairphi))
}


res2 <- fff2()

res2s0.00 <- fff2s() # 0% surrender penalty
res2s0.05 <- fff2s() # 5% surrender penalty
res2s0.15 <- fff2s() # 15% surrender penalty
res2s1.00 <- fff2s() # 100% surrender penalty


print(cbind(res2s1.00[[2]], res2s0.15[[2]], res2s0.05[[2]], res2s0.00[[2]]) * 10000, digits = 4)





gggs <- function(){
  gV <- seq(0.04, 0.06, by = 0.005)  
  
  fairphi <- vector(mode = "numeric", length = length(gV))
  
  for (i in 1 : length(gV))
  {
    
    phiV <- c(0, 0.5, seq(0.001, 0.01, by = 0.001), 0.05, 0.1)
    VV <- vector(mode = "numeric", length = length(phiV))
    
    for (j in 1 : length(phiV))
    {
      
      XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phiV[j], beta)
      
      FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phiV[j], beta = 0, gt = gV[i], tau, B)
      
      VV[j] <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 1)[[1]]
      
      print(c(gV[i], phiV[j], VV[j]))
      
      if (j == 1 & VV[1] < X0) {fairphi[i] <- -99; break()}
      if (j == 2 & VV[2] > X0) {fairphi[i] <- 99; break()}
      
    }
    if (fairphi[i] %in% c(-99, 99)) next()
    
    splineV <- function(fff)splinefun(x = phiV, y = VV, method = "hyman")(fff) - X0
    # splineV <- function(fff)approxfun(x = phiV, y = VV)(fff) - X0
    plot(splineV, from = min(phiV), to = max(phiV))
    fairphi[i] <- uniroot(f = splineV, interval = range(phiV))$root
    print(fairphi[i])
    
    XM <- simulate_XM(nsim = nsim, X0 = 1, M0 = X0, dS_S = dS_S, N = N, d = d, psi = psi, phi = fairphi[i], beta = beta)
    FD <- calculate_F(nsim = nsim, Xt = XM[[1]], Mt = XM[[2]], dS_S = dS_S, N = N, d = d, psi = psi, phi = fairphi[i], beta = 0, gt = gV[i], tau = tau, B = B)
    print(optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = 1)[[1]])
  }
  return(list(gV, fairphi))
}


gggs()





ppps <- function(){
  penaltyV <- c(seq(0, 0.1, by = 0.01), 0.2, 0.5, 1)
  VV <- vector(mode = "numeric", length = length(penaltyV))
  for (i in 1 : length(penaltyV))
  {
  XM <- simulate_XM(nsim, X0 = 1, M0 = X0, dS_S, N, d, psi, phi, beta)
  FD <- calculate_F(nsim, XM[[1]], XM[[2]], dS_S, N, d, psi, phi, beta = 0, gt = g, tau, B)
      
  VV[i] <- optimal_initiation_early_surrender(nsim = nsim, X = XM[[1]], M = XM[[2]], state_variables = NULL, F = FD, basis_functions = basis_functions_2v[[3]], n_basis_functions = 9, N = 60, d = 1, tau = tau, B = B, penalty = penaltyV[i])[[1]]
        
  print(c(penaltyV[i], VV[i]))
  }
  return(VV)
}


ppps()




