# 11. Update sig2_1
update_sig2_1 <- function(Y1, E1, F1, M1, X1, Z1, a_sig2_1, b_sig2_1, G, n) {
    sig2_1 = 1/rgamma(G, shape = (a_sig2_1 + n/2), scale = 1/(b_sig2_1 + rowSums((Y1 -
        M1 - E1 %*% X1 - F1 %*% Z1)^2)/2))
    return(sig2_1)
}
