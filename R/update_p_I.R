
update_p_I <- function(I, a_p_I, b_p_I) {
    # a_p_I and b_p_I are the hyperparameters
    G <- length(I)
    s <- sum(I)
    p_I <- rbeta(1, shape1 = s + a_p_I, shape2 = G - s + b_p_I)
    return(p_I)
}
