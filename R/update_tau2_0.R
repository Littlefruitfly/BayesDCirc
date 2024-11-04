# 14. Update tau2_0
update_tau2_0 <- function(I, mu_E, mu_F, a_tau2_0, b_tau2_0) {
    tau2_0 = 1/rgamma(1, shape = (a_tau2_0 + sum(I == 0)), rate = (b_tau2_0 + (sum((mu_E[I ==
        0])^2) + sum((mu_F[I == 0])^2))/2))
    return(tau2_0)
}
