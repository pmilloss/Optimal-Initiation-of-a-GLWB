simulate_tau_gompertz <- function(nsim, x0, b, m)
{
  b * log(rexp(nsim) * exp(-(x0 - m) / b) + 1)
}