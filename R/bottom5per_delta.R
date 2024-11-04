##' Get prior for delta
##'
##' We want to get prior for delta
##' @title Get prior for delta
##' @param mu_M mu_M from linear regression 
##' @return prior
##' @author Yutao Zhang
##' @export
##' @examples
##' mu_M = rnorm(1000)
##' bottom5per_delta(mu_M)

bottom5per_delta <- function(mu_M) {
    df <- data.frame(mu_M)
    df$square <- df$mu_M^2
    sorted_df <- df[order(df$square), ]
    sorted_df <- sorted_df[1:(0.05 * nrow(sorted_df)), ]
    a_delta2_0 <- (1/2) * nrow(sorted_df)
    b_delta2_0 <- (1/2) * sum(sorted_df$mu_M^2)
    prior <- data.frame(a_delta2_0, b_delta2_0)
    return(prior)
}
