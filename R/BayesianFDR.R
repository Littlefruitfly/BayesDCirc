##' BayesianFDR Function
##'
##' We want to compute the q values using this function
##' @title compute the q values
##' @param belief a vector, the ratio of sums of indicator for each gene to the iteration after burning 
##' @return q value
##' @author Yutao Zhang
##' @export
##' @examples
##' belief = runif(1000)
##' BayesianFDR(belief)

BayesianFDR <- function(belief) {
    pvalue <- 1 - belief
    pvalue_order <- order(pvalue)
    sortedP <- pvalue[pvalue_order]
    sortedQ <- cumsum(sortedP)/(1:length(sortedP))
    qvalue <- sortedQ[match(1:length(pvalue_order), pvalue_order)]
    return(qvalue)
}
