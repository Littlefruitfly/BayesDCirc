# 21. Update delta2_1
update_delta2_1 <- function(L, mu_M, a_delta2_1, b_delta2_1) {
    delta2_1 = 1/rgamma(1, shape = (a_delta2_1 + (1/2) * sum(L == 1)), rate = (b_delta2_1 +
        sum((mu_M[L == 1])^2)/2))
    return(delta2_1)
}
