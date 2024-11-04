# 6. Update coefficient F1
update_F1 <- function(Y1, Y2, X1, Z1, X2, Z2, M1, M2, sig2_F1, sig2_1, sig2_2, E1, E2,
    mu_F) {
    var_F = 1/(1/sig2_F1 + sum(Z1^2)/sig2_1 + sum(Z2^2)/sig2_2)
    mean_F = var_F * (((Y1 - M1 - E1 %*% X1) %*% t(Z1))/sig2_1 + ((Y2 - M2 - E2 %*%
        X2 - mu_F %*% Z2) %*% t(Z2))/sig2_2)
    F1 = rnorm(length(var_F), mean_F, sqrt(var_F))
    # use_F1 <- cbind(mean_F, var_F) F1 <- apply(use_F1, 1, function(x) rnorm(1,
    # x[1], sqrt(x[2])))
    return(F1)
}
