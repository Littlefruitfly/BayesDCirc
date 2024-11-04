##' gibbs sampling
##'
##' We want to use gibbs sampling to sample from posterior distributions for our parameter (Details)
##' @title gibbs sampling main function
##' @param Y1 expression values of group 1
##' @param Y2 expression values of group 2
##' @param t1 phenotype data of group 1
##' @param t2 phenotype data of group 2
##' @param n.iter number of iteration
##' @param hyperparameters hyperparameters we need
##' @param showIteration logic value, if need show iteration when processing gibbs sampling
##' @param priorInformation a data frame with three columns. The first is gene's name, the second is the indicator for differential circadian gene with pattern difference, and the third is the indicator for differential circadian gene with MESOR difference. Default is NULL.
##' @param priorFDR The FDR cut off that prior information was used to claim differential circadian genes. Default is NULL.
##' @return a list that contains posterior values
##' @author Yutao Zhang
##' @export
gibbs_main <- function(Y1, Y2, t1, t2, n.iter, hyperparameters, showIteration, priorInformation = NULL, priorFDR = NULL){
    # Y1 Y2 t1 t2 priorInformation should have the same row names
    t1 <- as.vector(t1)
    t2 <- as.vector(t2)
    t1 <- t(as.matrix(t1))
    t2 <- t(as.matrix(t2))
    X1 = cos((t1*pi)/12) # w = (2*pi)/Period
    Z1 = sin((t1*pi)/12)
    X2 = cos((t2*pi)/12)
    Z2 = sin((t2*pi)/12)
    
    G.record <- ceiling(n.iter/2) # ceiling: Round Up
    
    if(nrow(Y1) != nrow(Y2)){
        message("  Error: gene number is different between two group\n")
    }
    
    
    G <- nrow(Y1)
    G_d_I = 0
    G_d_L = 0
    n1 <- ncol(Y1)
    n2 <- ncol(Y2)
    
    if(!is.null(priorInformation)){
        if(!identical(rownames(Y1),as.character(priorInformation[,1]))){
            stop("prior Information doesn't match, check the row names of data")
        }
        if(is.null(priorFDR)){
            stop("  Please input priorFDR argument!\n")
        }
        G_d_I <- sum(priorInformation[,2])
        G_d_L <- sum(priorInformation[,3])
        #G_n <- G - sum(priorInformation[,2])
        index_I <- priorInformation[,2] == 1 # the I indicator index of diffcircadian genes in prior information
        index_L <- priorInformation[,3] == 1 # the L indicator index of diffcircadian genes in prior information 
    }
    
    #specify hyperparameters
    a_p_I <- hyperparameters[1]
    b_p_I <- hyperparameters[2]
    a_p_L <- hyperparameters[3]
    b_p_L <- hyperparameters[4]
    
    sig2_E1 <- hyperparameters[5]
    sig2_F1 <- hyperparameters[6]
    sig2_M1 <- hyperparameters[7]
    
    a_sig2_1 <- hyperparameters[8]
    b_sig2_1 <- hyperparameters[9]
    a_sig2_2 <- hyperparameters[10]
    b_sig2_2 <- hyperparameters[11]
    
    a_tau2_1 <- hyperparameters[12]
    b_tau2_1 <- hyperparameters[13]
    a_tau2_0 <- hyperparameters[14]
    b_tau2_0 <- hyperparameters[15]
    a_delta2_1 <- hyperparameters[16]
    b_delta2_1 <- hyperparameters[17]
    a_delta2_0 <- hyperparameters[18]
    b_delta2_0 <- hyperparameters[19]
    
    #set initial values
    # Setting up starting values
    I <- rep(0,G)
    L <- rep(0,G)
    E1 = rep(0.5,G)
    F1 = rep(0.5,G)
    M1 = rep(0.5,G)
    E2 = rep(0.5,G)
    F2 = rep(0.5,G)
    M2 = rep(0.5,G)
    mu_E = rep(0,G)
    mu_F = rep(0,G)
    mu_M = rep(0,G)
    sig2_1 <- rep(0.5, G)
    sig2_2 <- rep(0.5, G)
    tau2_0 <- 0.1
    tau2_1 <- 5
    delta2_0 <- 0.1
    delta2_1 <- 5
    
    #record samples
    p_I_record <- array(NA, dim = c(G.record))
    p_L_record <- array(NA, dim = c(G.record))
    I_record <- array(NA, dim = c(G, G.record))
    L_record <- array(NA, dim = c(G, G.record))
    E1_record <- array(NA, dim = c(G, G.record))
    F1_record <- array(NA, dim = c(G, G.record))
    M1_record <- array(NA, dim = c(G, G.record))
    mu_E_record <- array(NA, dim = c(G, G.record))
    mu_F_record <- array(NA, dim = c(G, G.record))
    mu_M_record <- array(NA, dim = c(G, G.record))
    E2_record <- array(NA, dim = c(G, G.record))
    F2_record <- array(NA, dim = c(G, G.record))
    M2_record <- array(NA, dim = c(G, G.record))
    sig2_1_record <- array(NA, dim = c(G, G.record))
    sig2_2_record <- array(NA, dim = c(G, G.record))
    tau2_1_record <- array(NA, dim = c(G.record))
    tau2_0_record <- array(NA, dim = c(G.record))
    delta2_1_record <- array(NA, dim = c(G.record))
    delta2_0_record <- array(NA, dim = c(G.record))
    
    #Gibbs sampler
    t1_sys <- Sys.time()
    message("  running the Gibbs sampler ...\n")
    for (i in 1:n.iter) {
        
        if(G_d_I != 0){
            p_I_d_with_prior <- update_p_I(I[index_I],G_d_I * (1 - priorFDR),G_d_I * priorFDR)
            p_I_n_with_prior <- update_p_I(I[!index_I],a_p_I,b_p_I)
            I[index_I] <- update_I(p_I_d_with_prior,mu_E[index_I],mu_F[index_I],tau2_1,tau2_0)
            I[!index_I] <- update_I(p_I_n_with_prior,mu_E[!index_I],mu_F[!index_I],tau2_1,tau2_0)
        }else{
            #Sample p_I
            p_I <- update_p_I(I,a_p_I,b_p_I)
            #Sample I
            I <- update_I(p_I,mu_E,mu_F,tau2_1,tau2_0)
        }
        
        if(G_d_L != 0){
            p_L_d_with_prior <- update_p_L(L[index_L],G_d_L * (1 - priorFDR), G_d_L * priorFDR)
            p_L_n_with_prior <- update_p_L(L[!index_L],a_p_L,b_p_L)
            L[index_L] <- update_L(p_L_d_with_prior,mu_M[index_L],delta2_1,delta2_0)
            L[!index_L] <- update_L(p_L_n_with_prior,mu_M[!index_L],delta2_1,delta2_0)
        }else{
            #Sample p_L
            p_L <- update_p_L(L,a_p_L,b_p_L)
            #Sample L
            L <- update_L(p_L,mu_M,delta2_1,delta2_0)
        }
        
        #Sample E1
        E1 <- update_E1(Y1, Y2, X1, Z1, X2, Z2, M1, M2, sig2_E1, sig2_1, sig2_2, F1, F2, mu_E)
        #Sample F1
        F1 <- update_F1(Y1, Y2, X1, Z1, X2, Z2, M1, M2, sig2_F1, sig2_1, sig2_2, E1, E2, mu_F)
        #Sample M1
        M1 <- update_M1(Y1, Y2, X1, Z1, X2, Z2, E1, E2, F1, F2, sig2_M1, sig2_1, sig2_2, mu_M,n1 = n1,n2 = n2)
        #Sample mu_E
        mu_E <- update_mu_E(Y1, Y2, X2, Z2, M1, M2, tau2_1, tau2_0, sig2_2, E1, F2, I)
        #Update E2
        E2 <- update_E2(E1,mu_E)
        #Sample mu_F
        mu_F <- update_mu_F(Y1, Y2, X2, Z2, M1, M2, tau2_1, tau2_0, sig2_2, E2, F1, I)
        #Update F2
        F2 <- update_F2(F1,mu_F)
        #Sample mu_M
        mu_M <- update_mu_M(Y1, Y2, X2, Z2, E2, F2, M1, delta2_1, delta2_0, sig2_2,L,n2=n2)
        #Update M2
        M2 <- update_M2(M1,mu_M)
        #Sample sig2
        sig2_1 <- update_sig2_1(Y1, E1, F1, M1, X1, Z1, a_sig2_1,b_sig2_1,G,n = n1)
        sig2_2 <- update_sig2_2(Y2, E2, F2, M2, X2, Z2, a_sig2_2,b_sig2_2,G,n = n2)
        #Sample tau2
        tau2_1 <- update_tau2_1(I,mu_E,mu_F,a_tau2_1,b_tau2_1)
        tau2_0 <- update_tau2_0(I,mu_E,mu_F,a_tau2_0,b_tau2_0)
        #Sample delta2
        delta2_1 <- update_delta2_1(L,mu_M,a_delta2_1,b_delta2_1)
        delta2_0 <- update_delta2_0(L,mu_M,a_delta2_0,b_delta2_0)  
        
        if(showIteration == TRUE){
            message(c("  Iteration ", i, "\n") )
        }
        
        ##store sampler
        if(i > n.iter - G.record){
            if(G_d_I != 0){
                p_I_record[i - (n.iter - G.record)] <- p_I_d_with_prior
            }else{
                p_I_record[i - (n.iter - G.record)] <- p_I
            }
            if(G_d_L != 0){
                p_L_record[i - (n.iter - G.record)] <- p_L_d_with_prior
            }else{
                p_L_record[i - (n.iter - G.record)] <- p_L
            }
            I_record[ ,i - (n.iter - G.record)] <- I
            L_record[ ,i - (n.iter - G.record)] <- L
            E1_record[ ,i - (n.iter - G.record)] <- E1
            F1_record[ ,i - (n.iter - G.record)] <- F1
            M1_record[ ,i - (n.iter - G.record)] <- M1
            mu_E_record[ ,i - (n.iter - G.record)] <- mu_E
            mu_F_record[ ,i - (n.iter - G.record)] <- mu_F
            mu_M_record[ ,i - (n.iter - G.record)] <- mu_M
            E2_record[ ,i - (n.iter - G.record)] <- E2
            F2_record[ ,i - (n.iter - G.record)] <- F2
            M2_record[ ,i - (n.iter - G.record)] <- M2
            sig2_1_record[ ,i - (n.iter - G.record)] <- sig2_1
            sig2_2_record[ ,i - (n.iter - G.record)] <- sig2_2
            tau2_1_record[i - (n.iter - G.record)] <- tau2_1
            tau2_0_record[i - (n.iter - G.record)] <- tau2_0
            delta2_1_record[i - (n.iter - G.record)] <- delta2_1
            delta2_0_record[i - (n.iter - G.record)] <- delta2_0
        }
    }
    t2_sys <- Sys.time()
    message(paste0("  The Gibbs sampler takes: ", round(difftime(t2_sys, t1_sys, units = "mins"), 3), " mins", "\n"))
    #posterior
    output <- list()
    
    ##posterior mean of other parameters
    p_I_post <- mean(p_I_record)
    p_L_post <- mean(p_L_record)
    E1_post <- rowMeans(E1_record)
    F1_post <- rowMeans(F1_record)
    M1_post <- rowMeans(M1_record)
    mu_E_post <- rowMeans(mu_E_record)
    mu_F_post <- rowMeans(mu_F_record)
    mu_M_post <- rowMeans(mu_M_record)
    E2_post <- rowMeans(E2_record)
    F2_post <- rowMeans(F2_record)
    M2_post <- rowMeans(M2_record)
    sig2_1_post <- rowMeans(sig2_1_record)
    sig2_2_post <- rowMeans(sig2_2_record)
    
    
    output[[1]] <- p_I_post
    output[[2]] <- p_L_post
    output[[3]] <- E1_post
    output[[4]] <- F1_post
    output[[5]] <- M1_post
    output[[6]] <- mu_E_post
    output[[7]] <- mu_F_post
    output[[8]] <- mu_M_post
    output[[9]] <- E2_post
    output[[10]] <- F2_post
    output[[11]] <- M2_post
    output[[12]] <- sig2_1_post
    output[[13]] <- sig2_2_post
    
    
    output[[14]] <- tau2_1_record
    output[[15]] <- tau2_0_record
    output[[16]] <- delta2_1_record
    output[[17]] <- delta2_0_record
    output[[18]] <- p_I_record
    output[[19]] <- p_L_record
    output[[20]] <- I_record
    output[[21]] <- L_record
    output[[22]] <- L_record|I_record
    
    output[[23]] <- round(difftime(t2_sys, t1_sys, units = "mins"), 3)
    
    output[[24]] <- E1_record
    output[[25]] <- F1_record
    output[[26]] <- M1_record
    output[[27]] <- E2_record
    output[[28]] <- F2_record
    output[[29]] <- M2_record
    
    names(output) <- c("p_I_post","p_L_post", "E1", "F1", "M1","mu_E", "mu_F", "mu_M", "E2", "F2","M2", "sig2_1_post", "sig2_2_post",
                       "tau2_1_record","tau2_0_record","delta2_1_record","delta2_0_record","p_I_record","p_L_record","I_record","L_record","df_circadian_indicators","computing_time",
                       "E1_record","F1_record","M1_record","E2_record","F2_record","M2_record")
    
    return(output)
}