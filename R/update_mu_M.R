# 19. Update difference mu_M
update_mu_M <- function(Y1, Y2, X2, Z2, E2, F2, M1, delta2_1, delta2_0, sig2_2, L, n2) {

    var_mu = 1/(1/delta2_1 + n2/sig2_2)
    mean_mu = var_mu * (rowSums(Y2 - E2 %*% X2 - F2 %*% Z2 - M1)/sig2_2)
    mu_M = rnorm(length(var_mu), mean_mu, sqrt(var_mu))

    var_mu = 1/(1/delta2_0 + n2/sig2_2)
    mean_mu = var_mu * (rowSums(Y2 - E2 %*% X2 - F2 %*% Z2 - M1)/sig2_2)
    mu_M_temp = rnorm(length(var_mu), mean_mu, sqrt(var_mu))

    mu_M[!L] = mu_M_temp[!L]
    return(mu_M)
}
