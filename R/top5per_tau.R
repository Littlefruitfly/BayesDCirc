##' Get prior for tau
##'
##' We want to get prior for tau
##' @title Get prior for tau
##' @param mu_E mu_E from linear regression 
##' @param mu_F mu_F from linear regression
##' @return prior
##' @author Yutao Zhang
##' @export
##' @examples
##' mu_E = rnorm(1000)
##' mu_F = rnorm(1000)
##' top5per_tau(mu_E,mu_F)
top5per_tau <- function(mu_E, mu_F) {
    df <- data.frame(mu_E, mu_F)
    df$square <- df$mu_E^2 + df$mu_F^2
    sorted_df <- df[order(-df$square), ]
    sorted_df <- sorted_df[1:(0.05 * nrow(sorted_df)), ]
    a_tau2_1 <- nrow(sorted_df)
    b_tau2_1 <- (1/2) * (sum(sorted_df$mu_E^2) + sum(sorted_df$mu_F^2))
    prior <- data.frame(a_tau2_1, b_tau2_1)
    return(prior)
}
