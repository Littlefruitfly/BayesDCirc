# 13. Update tau2_1
update_tau2_1 <- function(I, mu_E, mu_F, a_tau2_1, b_tau2_1) {
    tau2_1 = 1/rgamma(1, shape = (a_tau2_1 + sum(I == 1)), rate = (b_tau2_1 + (sum((mu_E[I ==
        1])^2) + sum((mu_F[I == 1])^2))/2))
    return(tau2_1)
}
