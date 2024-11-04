# 18. Update mesor M1
update_M1 <- function(Y1, Y2, X1, Z1, X2, Z2, E1, E2, F1, F2, sig2_M1, sig2_1, sig2_2,
    mu_M, n1, n2) {
    var_M = 1/(1/sig2_M1 + n1/sig2_1 + n2/sig2_2)
    mean_M = var_M * (rowSums(Y1 - E1 %*% X1 - F1 %*% Z1)/sig2_1 + rowSums(Y2 - E2 %*%
        X2 - F2 %*% Z2 - mu_M)/sig2_2)
    M1 = rnorm(length(var_M), mean_M, sqrt(var_M))
    return(M1)
}
