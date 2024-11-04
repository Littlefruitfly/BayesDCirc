# 12. Update sig2_2
update_sig2_2 <- function(Y2, E2, F2, M2, X2, Z2, a_sig2_2, b_sig2_2, G, n) {
    sig2_2 = 1/rgamma(G, shape = (a_sig2_2 + n/2), scale = 1/(b_sig2_2 + rowSums((Y2 -
        M2 - E2 %*% X2 - F2 %*% Z2)^2)/2))
    return(sig2_2)
}
