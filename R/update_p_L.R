# 16. Update the differentical MESOR proportion
update_p_L <- function(L, a_p_L, b_p_L) {
    # a_p and b_p are the hyperparameters for
    G <- length(L)
    s <- sum(L)
    p_L <- rbeta(1, shape1 = s + a_p_L, shape2 = G - s + b_p_L)
    return(p_L)
}
