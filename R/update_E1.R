# 5. Update coefficient E1
update_E1 <- function(Y1, Y2, X1, Z1, X2, Z2, M1, M2, sig2_E1, sig2_1, sig2_2, F1, F2,
    mu_E) {
    var_E = 1/(1/sig2_E1 + sum(X1^2)/sig2_1 + sum(X2^2)/sig2_2)
    mean_E = var_E * (((Y1 - M1 - F1 %*% Z1) %*% t(X1))/sig2_1 + ((Y2 - M2 - F2 %*%
        Z2 - mu_E %*% X2) %*% t(X2))/sig2_2)
    E1 = rnorm(length(var_E), mean_E, sqrt(var_E))
    # use_E1 <- cbind(mean_E, var_E) E1 <- apply(use_E1, 1, function(x) rnorm(1,
    # x[1], sqrt(x[2])))
    return(E1)
}
