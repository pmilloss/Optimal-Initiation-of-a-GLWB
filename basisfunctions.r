# polynomials basis functions for two variables, polynomial degree 1:4
basis_functions_2v <- list(function(x)x,
                           function(x)c(x, x ^ 2, x[1] * x[2]),
                           function(x)c(x, x ^ 2, x[1] * x[2], x ^ 3, x[1] ^ 2 * x[2], x[1] * x[2] ^ 2),
                           function(x, y)c(x, x ^ 2, x[1] * x[2], x ^ 3, x[1] ^ 2 * x[2], x[1] * x[2] ^ 2, x ^ 4, x[1] * x[2] ^ 3, x[1] ^ 2 * x[2] ^ 2, x[1] ^ 3 * x[2]))# 2, 5, 9, 14

# polynomials basis functions for three variables, polynomial degree 1:4
basis_functions_3v <- list(function(x)x,
                        function(x)c(x, x ^ 2, x[1] * x[2], x[1] * x[3], x[2] * x[3]),
                        function(x)c(x, x ^ 2, x[1] * x[2], x[1] * x[3], x[2] * x[3], x ^ 3, x[1] ^ 2 * x[2], x[1] * x[2] ^ 2, x[1] * x[3] ^ 2, x[1] ^ 2 * x[3], x[2] * x[3] ^ 2, x[2] ^ 2 * x[3], x[1] * x[2] * x[3]),
                        function(x)c(x, x ^ 2, x[1] * x[2], x[1] * x[3], x[2] * x[3], x ^ 3, x[1] ^ 2 * x[2], x[1] * x[2] ^ 2, x[1] * x[3] ^ 2, x[1] ^ 2 * x[3], x[2] * x[3] ^ 2, x[2] ^ 2 * x[3], x[1] * x[2] * x[3], x ^ 4, x[1] * x[2] ^ 3, x[1] ^ 2 * x[2] ^ 2, x[1] ^ 3 * x[2], x[1] * x[3] ^ 3, x[1] ^ 2 * x[3] ^ 2, x[1] ^ 3 * x[3], x[2] * x[3] ^ 3, x[2] ^ 2 * x[3] ^ 2, x[2] ^ 3 * x[3], x[1] ^ 2 * x[2] * x[3], x[1] * x[2] ^ 2 * x[3], x[1] * x[2] * x[3] ^ 2)
)# 3, 9, 19, 34

# polynomials basis functions for four variables, polynomial degree 1:4
basis_functions_4v <- list(function(x)x,
                     function(x)c(x, x ^ 2, x[1] * x[2], x[1] * x[3], x[1] * x[4], x[2] * x[4], x[2] * x[4], x[3] * x[4]),
                     function(x)c(x, x ^ 2, x[1] * x[2], x[1] * x[3], x[1] * x[4], x[2] * x[3], x[2] * x[4], x[3] * x[4], x ^ 3, x[1] ^ 2 * x[2], x[1] * x[2] ^ 2, x[1] * x[3] ^ 2, x[1] ^ 2 * x[3], x[1] * x[4] ^ 2, x[1] ^ 2 * x[4], x[2] * x[3] ^ 2, x[2] ^ 2 * x[3], x[2] * x[4] ^ 2, x[2] ^ 2 * x[4], x[3] * x[4] ^ 2, x[3] ^ 2 * x[4],
x[1] * x[2] * x[3], x[1] * x[2] * x[4], x[1] * x[3] * x[4], x[2] * x[3] * x[4]),
                     function(x)c(x, x ^ 2, x[1] * x[2], x[1] * x[3], x[1] * x[4], x[2] * x[3], x[2] * x[4], x[3] * x[4], x ^ 3, x[1] ^ 2 * x[2], x[1] * x[2] ^ 2, x[1] * x[3] ^ 2, x[1] ^ 2 * x[3], x[1] * x[4] ^ 2, x[1] ^ 2 * x[4], x[2] * x[3] ^ 2, x[2] ^ 2 * x[3], x[2] * x[4] ^ 2, x[2] ^ 2 * x[4], x[3] * x[4] ^ 2, x[3] ^ 2 * x[4], x[1] * x[2] * x[3], x[1] * x[2] * x[4], x[1] * x[3] * x[4], x[2] * x[3] * x[4], x ^ 4, x[1] * x[2] ^ 3, x[1] ^ 2 * x[2] ^ 2, x[1] ^ 3 * x[2], x[1] * x[3] ^ 3, x[1] ^ 2 * x[3] ^ 2, x[1] ^ 3 * x[3], x[1] * x[4] ^ 3, x[1] ^ 2 * x[4] ^ 2, x[1] ^ 3 * x[4], x[2] * x[3] ^ 3, x[2] ^ 2 * x[3] ^ 2, x[2] ^ 3 * x[3], x[2] * x[4] ^ 3, x[2] ^ 2 * x[4] ^ 2, x[2] ^ 3 * x[4], x[3] * x[4] ^ 3, x[3] ^ 2 * x[4] ^ 2, x[3] ^ 3 * x[4], x[1] ^ 2 * x[2] * x[3], x[1] * x[2] ^ 2 * x[3], x[1] * x[2] * x[3] ^ 2, x[1] ^ 2 * x[2] * x[4], x[1] * x[2] ^ 2 * x[4], x[1] * x[2] * x[4] ^ 2, x[1] ^ 2 * x[3] * x[4], x[1] * x[3] ^ 2 * x[4], x[1] * x[3] * x[4] ^ 2, x[2] ^ 2 * x[3] * x[4], x[2] * x[3] ^ 2 * x[4], x[2] * x[3] * x[4] ^ 2, x[1] * x[2] * x[3] * x[4])
)# 4, 14, 34, 69