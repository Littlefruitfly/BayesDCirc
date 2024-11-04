
update_L <- function(p_L, mu_M, delta2_1, delta2_0) {
    # output is a vector with G components
    prob1 <- log(p_L) + dnorm(mu_M, 0, sqrt(delta2_1), log = T)
    prob2 <- log(1 - p_L) + dnorm(mu_M, 0, sqrt(delta2_0), log = T)
    L <- rbinom(n = length(prob1), size = 1, prob = 1/(1 + exp(prob2 - prob1)))
    return(L)
}
