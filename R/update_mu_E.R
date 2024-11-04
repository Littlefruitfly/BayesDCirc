# 7. Update difference mu_E
update_mu_E <- function(Y1, Y2, X2, Z2, M1, M2, tau2_1, tau2_0, sig2_2, E1, F2, I) {

    var_mu = 1/(1/tau2_1 + sum(X2^2)/sig2_2)
    mean_mu = var_mu * (((Y2 - M2 - F2 %*% Z2 - E1 %*% X2) %*% t(X2))/sig2_2)
    mu_E = rnorm(length(var_mu), mean_mu, sqrt(var_mu))

    var_mu = 1/(1/tau2_0 + sum(X2^2)/sig2_2)
    mean_mu = var_mu * (((Y2 - M2 - F2 %*% Z2 - E1 %*% X2) %*% t(X2))/sig2_2)
    mu_E_temp = rnorm(length(var_mu), mean_mu, sqrt(var_mu))

    # mu_E[I] = mu_E_temp1[I]
    mu_E[!I] = mu_E_temp[!I]
    return(mu_E)
}
