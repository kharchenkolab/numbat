#'@import Rcpp
NULL

############ Allele HMMs ############

#' Get an allele HMM
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta numeric Haplotype imbalance
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @return HMM object
#' @keywords internal
get_allele_hmm = function(pAD, DP, p_s, theta, gamma = 20) {

    states = c("theta_up", "theta_down")

    N = length(p_s)

    Pi = sapply(p_s, function(p_s) {c(1 - p_s, p_s, p_s, 1 - p_s)}) %>% 
        array(dim = c(2, 2, N))

    if (length(theta) == 1) {
        theta = rep(theta, N)
    }

    prior = c(0.5, 0.5)
    alpha_up = (0.5 + theta) * gamma
    beta_up = (0.5 - theta) * gamma
    alpha_down = beta_up
    beta_down = alpha_up
    
    hmm = list(
        x = pAD, 
        logPi = log(Pi),
        delta = prior, 
        alpha = matrix(c(alpha_up, alpha_down), ncol = 2), 
        beta = matrix(c(beta_up, beta_down), ncol = 2),
        d = DP,
        N = N,
        M = 2,
        K = 1,
        states = states
    )

    return(hmm)
}


#' Calculate allele likelihoods
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta numeric Haplotype imbalance
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @keywords internal                         
calc_allele_lik = function (pAD, DP, p_s, theta, gamma = 20) {
    hmm = get_allele_hmm(pAD, DP, p_s, theta, gamma)
    LL = likelihood_allele(hmm)
    return(LL)
}

############ Clonal deletion HMM ############

#' Viterbi for clonal LOH detection
#' @param hmm HMM object; expect variables x (SNP count), snp_sig (snp rate standard deviation), 
#' pm (snp density for ref and loh states), pn (gene lengths),
#' d (total expression depth), y (expression count), lambda_star (reference expression rate), 
#' mu (global expression mean), sig (global expression standard deviation),
#' Pi (transition prob matrix), delta (prior for each state),
#' phi (expression fold change for each state)
#' @keywords internal
viterbi_loh <- function (hmm, ...){

    n <- length(hmm$x)
    m <- nrow(hmm$Pi[,,1])
    nu <- matrix(NA, nrow = n, ncol = m)
    mu <- matrix(NA, nrow = n, ncol = m + 1)
    z <- rep(NA, n)
    
    nu[1, ] = log(hmm$delta)
    logPi <- log(hmm$Pi)

    for (i in 1:n) {

        if (i > 1) {
            matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
            nu[i, ] = apply(matrixnu + logPi[,,i], 2, max)
        }
        
        nu[i, ] = nu[i, ] + dnbinom(x = hmm$x[i], mu = hmm$pm * hmm$pn[i], size = hmm$snp_sig, log = TRUE)

        nu[i, ] = nu[i, ] + dpoilog(
            x = rep(hmm$y[i], m),
            sig = rep(hmm$sig, m),
            mu = hmm$mu + log(hmm$phi * hmm$d * hmm$lambda_star[i]),
            log = TRUE
        )
    }
             
    z[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1)) z[i] <- which.max(logPi[,,i+1][, z[i+1]] + nu[i, ])

    LL = max(nu[n, ])
        
    return(hmm$states[z])
}