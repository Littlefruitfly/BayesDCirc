update_I <- function(p_I, mu_E, mu_F, tau2_1, tau2_0) {
    # update I, I = 1 if the expression level of gene g is differential between
    # two groups output is a vector with G components
    prob1 <- log(p_I) + dnorm(mu_E, 0, sqrt(tau2_1), log = T) + dnorm(mu_F, 0, sqrt(tau2_1),
        log = T)
    prob2 <- log(1 - p_I) + dnorm(mu_E, 0, sqrt(tau2_0), log = T) + dnorm(mu_F, 0, sqrt(tau2_0),
        log = T)
    I <- rbinom(n = length(prob1), size = 1, prob = 1/(1 + exp(prob2 - prob1)))
    return(I)
}
