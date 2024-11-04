# 8. Update difference mu_F
update_mu_F <- function(Y1, Y2, X2, Z2, M1, M2, tau2_1, tau2_0, sig2_2, E2, F1, I) {

    var_mu = 1/(1/tau2_1 + sum(Z2^2)/sig2_2)
    mean_mu = var_mu * (((Y2 - M2 - E2 %*% X2 - F1 %*% Z2) %*% t(Z2))/sig2_2)
    mu_F = rnorm(length(var_mu), mean_mu, sqrt(var_mu))

    var_mu = 1/(1/tau2_0 + sum(Z2^2)/sig2_2)
    mean_mu = var_mu * (((Y2 - M2 - E2 %*% X2 - F1 %*% Z2) %*% t(Z2))/sig2_2)
    mu_F_temp = rnorm(length(var_mu), mean_mu, sqrt(var_mu))

    mu_F[!I] = mu_F_temp[!I]
    return(mu_F)
}
