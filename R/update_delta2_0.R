# 22. Update delta2_0
update_delta2_0 <- function(L, mu_M, a_delta2_0, b_delta2_0) {
    delta2_0 = 1/rgamma(1, shape = (a_delta2_0 + (1/2) * sum(L == 0)), rate = (b_delta2_0 +
        sum((mu_M[L == 0])^2)/2))
    return(delta2_0)
}
