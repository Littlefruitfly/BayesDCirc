##' get_hyperparameters Function
##'
##' We want to get hyperparameters using this function
##' @title gibbs sampling main function
##' @param Y1 expression values of group 1
##' @param Y2 expression values of group 2
##' @param t1 phenotype data of group 1
##' @param t2 phenotype data of group 2
##' @param use_default logic value, if the default/non-informative values for the prior are used. It does not implement a spike-and-slab prior. Default is TRUE.
##' @param personalized_para numeric vector, required if use_default is FALSE. Personalized parameters can be input here. Default is NULL.
##' @param Period numeric value. The time period for which the data needs to be estimated. Default is 24.
##' @return a numeric vector that contains hyperparameters that can be used for `BayesDCirc`.
##' @author Yutao Zhang
##' @export
##' @examples
##' hyperparameters <- get_hyperparameters(Y1,Y2,t1,t2, use_default = T)

get_hyperparameters <- function(Y1,Y2,t1,t2, use_default = T, personalized_para = NULL, Period = 24) {
  stopifnot(nrow(Y1) == nrow(Y2))
  G <- nrow(Y1)
  Y1 <- as.matrix(Y1)
  Y2 <- as.matrix(Y2)
  w = (2*pi)/Period
  X1 = cos(t1*w) 
  Z1 = sin(t1*w)
  X2 = cos(t2*w)
  Z2 = sin(t2*w)
  
  #1. fit a linear model to check how large the E, F and M are
  fitted_values1 <- matrix(NA,G,3)
  colnames(fitted_values1) <- c('M1','E1','F1')
  for (i in 1:G) {
    mod <- lm(Y1[i,] ~ as.vector(X1) + as.vector(Z1))
    fitted_values1[i,] <- mod$coefficients
  }
  
  fitted_values2 <- matrix(NA,G,3)
  colnames(fitted_values2) <- c('M2','E2','F2')
  for (i in 1:G) {
    mod <- lm(Y2[i,] ~ as.vector(X2) + as.vector(Z2))
    fitted_values2[i,] <- mod$coefficients
  }
  
  difference <- fitted_values2 - fitted_values1
  colnames(difference) <- c("diff_M","mu_E","mu_F")
  difference <- as.data.frame(difference)
  
  #2. get prior for tau2_1 and tau2_0
  prior <- top5per_tau(difference$mu_E,difference$mu_F)
  a_tau2_1 <- prior$a_tau2_1
  b_tau2_1 <- prior$b_tau2_1
  a_tau2_1 
  b_tau2_1
  
  prior <- bottom5per_tau(difference$mu_E,difference$mu_F)
  a_tau2_0 <- prior$a_tau2_0
  b_tau2_0 <- prior$b_tau2_0
  a_tau2_0
  b_tau2_0
  
  #3. get prior for delta2_1 and delta2_0
  prior <- top5per_delta(difference$diff_M)
  a_delta2_1 <- prior$a_delta2_1
  b_delta2_1 <- prior$b_delta2_1
  a_delta2_1
  b_delta2_1
  
  prior <- bottom5per_delta(difference$diff_M)
  a_delta2_0 <- prior$a_delta2_0
  b_delta2_0 <- prior$b_delta2_0
  a_delta2_0
  b_delta2_0
  
  
  #4. set hyperparameters
  ifelse(use_default, hyperparameters <- c(1,1,1,1, 10,10,10, 0.001,0.001,0.001,0.001,
                                           a_tau2_1,b_tau2_1,a_tau2_0,b_tau2_0,a_delta2_1,b_delta2_1,a_delta2_0,b_delta2_0),
         hyperparameters <- personalized_para)
  return(hyperparameters)
}
