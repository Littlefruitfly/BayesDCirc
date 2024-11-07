##' BayesianFDR Function
##'
##' We want to compute FDR using this function
##' @title compute the Bayesian FDR
##' @param output_BayesDCirc a list containing posterior values, typically obtained from `BayesDCirc`.
##' @param indicator a character string, either "I" or "L", indicating which specific indicator to check.
##' @return a data frame containing FDR of all gene candidates.
##' @author Yutao Zhang
##' @export
##' @examples
##' BayesianFDR(output, indicator = "I")

BayesianFDR <- function(output_BayesDCirc, indicator) {
    if(indicator == "I"){belief <- rowSums(output_BayesDCirc$I_record)/dim(output_BayesDCirc$I_record)[2]}
    if(indicator == "L"){belief <- rowSums(output_BayesDCirc$L_record)/dim(output_BayesDCirc$I_record)[2]}
    pvalue <- 1 - belief
    pvalue_order <- order(pvalue)
    sortedP <- pvalue[pvalue_order]
    sortedQ <- cumsum(sortedP)/(1:length(sortedP))
    qvalue <- sortedQ[match(1:length(pvalue_order), pvalue_order)]
    result_diff <- data.frame(gene_candidates = output_BayesDCirc$gene_candidates,
                              FDR = qvalue,
                              indicator = indicator)
    return(result_diff)
}
