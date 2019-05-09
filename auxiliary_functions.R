# extract submatrix from matrix x containing rows i and columns j - to be used with lapply
extract_row_col <- function(x, row, col)x[row, col]

trap <- function(f, d) 0.5 * d * sum(f[-length(f)] + f[-1])

CIR1 <- function(CIR0, k, theta, sigma, d, dW)pmax(k * theta * d + (1 - k * d) * CIR0 + sigma * sqrt(CIR0) * dW, 0)
