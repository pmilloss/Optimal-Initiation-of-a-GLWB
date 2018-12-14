CIR1 <- function(CIR0, k, theta, sigma, d, dW)pmax(k*theta*d+(1-k*d)*CIR0+sigma*sqrt(CIR0)*dW, 0)

dS1.S0 <- function(r0, V0, d, dWV, dWS, rhoSV)pmax(r0*d+sqrt(V0)*(rhoSV*dWV+sqrt(1-rhoSV^2)*dWS), -0.999)

S1 <- function(S0, r0, V0, d, dWV, dWS, rhoSV)pmax(S0*(1+r0*d+sqrt(V0)*(rhoSV*dWV+sqrt(1-rhoSV^2)*dWS)), 10E-5)

A1 <- function(A0, r0, V0, d, dWV, dWA, rhoAV, phi)pmax(A0*(1+(r0-phi)*d+sqrt(V0)*(rhoAV*dWV+sqrt(1-rhoAV^2)*dWA)), 0)