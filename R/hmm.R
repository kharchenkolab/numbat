#' Beta-binomial distribution density function
#' A distribution is beta-binomial if p, the probability of success, 
#' in a binomial distribution has a beta distribution with shape 
#' parameters alpha > 0 and beta > 0
#' For more details, see extraDistr::dbbinom
#'
#' @param x vector of quantiles
#' @param size number of trials (zero or more)
#' @param alpha numeric (default=1)
#' @param beta numeric (default=1)
#' @param log boolean (default=FALSE)
#' @return density values returned as numeric vector
#' @keywords internal
dbbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
    cppdbbinom(x, size, alpha, beta, log[1L])
}

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

#' Viterbi algorithm for allele HMM
#' @param hmm HMM object; expect variables x (allele depth), d (total depth),
#' logPi (log transition prob matrix), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), 
#' states (states), p_s (phase switch probs)
#' @keywords internal
viterbi_allele <- function(hmm) {

    N <- hmm$N
    M <- hmm$M
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)

    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })

    z = viterbi_compute(log(hmm$delta), logprob, hmm$logPi, N, M, nu, z)
        
    return(z)
}

#' Forward-backward algorithm for allele HMM
#' @param hmm HMM object; expect variables x (allele depth), d (total depth),
#' logPi (log transition prob matrix), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), 
#' states (states), p_s (phase switch probs)
#' @keywords internal
forward_back_allele = function(hmm) {

    # case of one-data point
    if (hmm$N == 1) {
        return(NA)
    }

    N <- hmm$N
    M <- hmm$M
        
    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })
        
    logphi <- log(hmm$delta)
        
    marginals = forward_backward_compute(logphi, logprob, hmm$logPi, N, M)

    colnames(marginals) = hmm$states

    return(marginals)
}

# Only compute total log likelihood
#' @param hmm HMM object; expect variables x (allele depth), d (total depth),
#' logPi (log transition prob matrix), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), 
#' states (states), p_s (phase switch probs)
#' @keywords internal
likelihood_allele = function(hmm) {
        
    N <- hmm$N
    M <- hmm$M
        
    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[,m], beta = hmm$beta[,m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })
        
    logphi <- log(hmm$delta)
        
    LL <- likelihood_compute(logphi, logprob, hmm$logPi, N, M)

    return(LL)
}

#' allele-only HMM
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param t numeric Transition probability between copy number states
#' @param theta_min numeric Minimum haplotype frequency deviation threshold
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @return character vector Decoded states
#' @keywords internal
run_allele_hmm = function(pAD, DP, p_s, t = 1e-5, theta_min = 0.08, gamma = 20, prior = NULL) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('More than one gamma parameter')
    }

    # states
    states = c("neu", "theta_1_up", "theta_1_down", "theta_2_up", "theta_2_down")

    N = length(pAD)
    M = 5

    # transition matrices
    calc_trans_mat = function(p_s, t) {
        matrix(
            c(1-t, t/4, t/4, t/4, t/4, 
             t/2, (1-t)*(1-p_s), (1-t)*p_s, t/4, t/4,
             t/2, (1-t)*p_s, (1-t)*(1-p_s), t/4, t/4,
             t/2, t/4, t/4, (1-t)*(1-p_s), (1-t)*p_s,
             t/2, t/4, t/4, (1-t)*p_s, (1-t)*(1-p_s)),
            ncol = 5,
            byrow = TRUE
        )
    }

    Pi = sapply(p_s, function(p_s) {calc_trans_mat(p_s, t)}) %>% 
        array(dim = c(M, M, N))

    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/5, 5)
    }

    theta_1 = theta_min
    theta_2 = 0.4
    alphas = gamma * c(0.5, 0.5 + theta_1, 0.5 - theta_1, 0.5 + theta_2, 0.5 - theta_2)
    betas = gamma * c(0.5, 0.5 - theta_1, 0.5 + theta_1, 0.5 - theta_2, 0.5 + theta_2)
            
    hmm = list(
        x = pAD, 
        logPi = log(Pi), 
        delta = prior, 
        alpha = matrix(rep(alphas, N), ncol = M, byrow = TRUE),
        beta = matrix(rep(betas, N), ncol = M, byrow = TRUE),
        d = DP,
        N = N,
        M = M,
        states = states
    )
    
    mpc = states[viterbi_allele(hmm)]
    
    return(mpc)
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

############ Joint HMM ############


#' Run joint HMM on a pseudobulk profile
#' @param pAD integer vector Paternal allele counts
#' @param DP integer vector Total alelle counts
#' @param p_s numeric vector Phase switch probabilities
#' @param theta_min numeric Minimum haplotype imbalance threshold
#' @param gamma numeric Overdispersion in the allele-specific expression
#' @param Y_obs numeric vector Observed gene counts
#' @param lambda_ref numeric vector Reference expression rates
#' @param d_total integer Total library size for expression counts
#' @param phi_del numeric Expected fold change for deletion
#' @param phi_amp numeric Expected fold change for amplification
#' @param phi_bamp numeric Expected fold change for balanced amplification
#' @param phi_bdel numeric Expected fold change for balanced deletion
#' @param mu numeric Global expression bias
#' @param sig numeric Global expression variance
#' @param t numeric Transition probability between copy number states
#' @param exp_only logical Whether to only use expression data
#' @param allele_only logical Whether to only use allele data
#' @keywords internal
run_joint_hmm = function(
    pAD, DP, p_s, Y_obs = 0, lambda_ref = 0, d_total = 0, theta_min = 0.08, theta_neu = 0,
    bal_cnv = TRUE, phi_del = 2^(-0.25), phi_amp = 2^(0.25), phi_bamp = phi_amp, phi_bdel = phi_del, 
    alpha = 1, beta = 1, 
    mu = 0, sig = 1,
    t = 1e-5, gamma = 18,
    prior = NULL, exp_only = FALSE, allele_only = FALSE,
    classify_allele = FALSE, phasing = TRUE, debug = FALSE
) {

    # states
    states = c(
        "1" = "neu", "2" = "del_1_up", "3" = "del_1_down", "4" = "del_2_up", "5" = "del_2_down",
        "6" = "loh_1_up", "7" = "loh_1_down", "8" = "loh_2_up", "9" = "loh_2_down", 
        "10" = "amp_1_up", "11" = "amp_1_down", "12" = "amp_2_up", "13" = "amp_2_down", 
        "14" = "bamp", "15" = "bdel"
    )

    states_cn = str_remove(states, '_up|_down')
    states_phase = str_extract(states, 'up|down')

    # relative abundance of states
    w = c('neu' = 1, 'del_1' = 1, 'del_2' = 1e-10, 'loh_1' = 1, 'loh_2' = 1e-10, 'amp_1' = 1, 'amp_2' = 1e-10, 'bamp' = 1e-4, 'bdel' = 1e-10)
        
    # intitial probabilities
    if (is.null(prior)) {
        # encourage CNV from telomeres
        prior = sapply(1:length(states), function(to){
                get_trans_probs(
                    t = min(t * 100, 1), p_s = 0, w,
                    cn_from = 'neu', phase_from = NA,
                    cn_to = states_cn[to], phase_to = states_phase[to])
            })
    }

    # to do: renormalize the probabilities after deleting states
    states_index = 1:length(states)

    if (!bal_cnv) {
        states_index = 1:13
    }
        
    if (exp_only) {
        pAD = rep(NA, length(pAD))
        p_s = rep(0, length(p_s))
    }
    
    if (allele_only) {
        states_index = c(1, 6:9)

        Y_obs = rep(NA, length(Y_obs))
    }

    if (!phasing) {
        states_index = c(1, 6)
        
        p_s = ifelse(is.na(pAD), p_s, 0)
        pAD = ifelse(pAD > (DP - pAD), pAD, DP - pAD)
        theta_neu = 0.1
        theta_min = 0.45
    }

    if (classify_allele) {
        states_index = c(6,7)
    }
    
    # transition matrices
    As = calc_trans_mat(t, p_s, w, states_cn, states_phase)

    theta_u_1 = 0.5 + theta_min
    theta_d_1 = 0.5 - theta_min

    theta_u_2 = 0.9
    theta_d_2 = 0.1

    theta_u_neu = 0.5 + theta_neu
    theta_d_neu = 0.5 - theta_neu

    # parameters for each state
    alpha_states = gamma * c(theta_u_neu, rep(c(theta_u_1, theta_d_1, theta_u_2, theta_d_2), 3), theta_u_neu, theta_u_neu)
    beta_states = gamma * c(theta_d_neu, rep(c(theta_d_1, theta_u_1, theta_d_2, theta_u_2), 3), theta_d_neu, theta_d_neu)
    phi_states = c(1, rep(phi_del, 2), rep(0.5, 2), rep(1, 4), rep(phi_amp, 2), rep(2.5, 2), phi_bamp, phi_bdel)
    
    # subset for relevant states
    prior = prior[states_index]
    As = As[states_index, states_index,]
    alpha_states = alpha_states[states_index]
    beta_states = beta_states[states_index]
    phi_states = phi_states[states_index]
    states = states[states_index] %>% setNames(1:length(.))
                
    N = length(Y_obs)

    if (length(mu) == 1) {
        mu = rep(mu, N)
        sig = rep(sig, N)
    }

    if (length(d_total) == 1) {
        d_total = rep(d_total, N)
    }

    hmm = list(
        x = matrix(pAD), 
        d = matrix(DP),
        y = matrix(Y_obs),
        l = matrix(d_total),
        lambda = matrix(lambda_ref),
        mu = matrix(mu),
        sig = matrix(sig),
        logPi = log(As), 
        phi = matrix(phi_states),
        delta = matrix(prior), 
        alpha = matrix(alpha_states),
        beta = matrix(beta_states),
        states = matrix(states),
        p_s = p_s
    )

    MPC = states[viterbi_joint(hmm)]
        
    return(MPC)
}


#' Calculate the transition matrix for joint HMM
#' @param t numeric; CNV state transition probability
#' @param p_s numeric; phase switch probability
#' @param w numeric; relative abundance of states
#' @param states_cn character; CNV states
#' @param states_phase character; haplotype phase states
#' @return array; transition matrix
#' @keywords internal
calc_trans_mat = function(t, p_s, w, states_cn, states_phase) {

    sapply(1:length(states_cn), function(from) {
        sapply(1:length(states_cn), function(to) {
            get_trans_probs(t, p_s, w, states_cn[from], states_phase[from], states_cn[to], states_phase[to])
        }) %>% t
    }) %>% t %>%
    array(dim = c(length(states_cn), length(states_cn), length(p_s)))

}


#' Helper function to calculate transition porbabilities
#' cn/phase are sclars, only p_s is vectorized
#' @param t numeric; CNV state transition probability
#' @param p_s numeric; phase switch probability
#' @param w numeric; relative abundance of states
#' @param cn_from character; originating CNV state
#' @param phase_from character; originating haplotype phase state
#' @param cn_to character; destination CNV state
#' @param phase_to character; destination haplotype phase state
#' @return numeric; transition probability
#' @keywords internal
get_trans_probs = function(t, p_s, w, cn_from, phase_from, cn_to, phase_to) {

    if (cn_from == 'neu' & cn_to == 'neu') {
        p_s = rep(0.5, length(p_s))
    }

    if (cn_from == cn_to) {
        if (is.na(phase_from) & is.na(phase_to)) {
            p = 1-t
            p = rep(p, length(p_s))
        } else if (phase_from == phase_to) {
            p = (1-t) * (1-p_s)
        } else {
            p = (1-t) * p_s
        }
    } else {
        p = t * w[[cn_to]]/sum(w[names(w)!=cn_from])
        if (!is.na(phase_to)) {
            p = p/2
        }
        p = rep(p, length(p_s))
    }
    
    return(p)
}

#' Generalized viterbi algorithm for joint HMM
#' @param hmm HMM object; expect variables x (allele depth), d (total depth), y (expression count), l (cell total library size),
#' lambda (reference expression rate), mu (global expression mean), sig (global expression standard deviation), 
#' logPi (log transition prob matrix), phi (expression fold change for each state), delta (prior for each state), 
#' alpha (alpha for each state), beta (beta for each state), states (states), p_s (phase switch probs)
#' @keywords internal
viterbi_joint <- function(hmm) {

    N <- nrow(hmm$x)
    M <- nrow(hmm$logPi[,,1])
    K <- ncol(hmm$x)
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)
    
    logprob = sapply(1:M, function(m) {

        sapply(
            1:K,
            function(k) {

                l_x = dbbinom(x = hmm$x[,k], size = hmm$d[,k], alpha = hmm$alpha[m,k], beta = hmm$beta[m,k], log = TRUE)

                l_x[is.na(l_x)] = 0

                if (!is.null(hmm$y)) {
                    l_y = rep(0,N)
                    valid = !is.na(hmm$y[,k])

                    l_y[valid] = dpoilog(
                            x = hmm$y[valid,k],
                            sig = hmm$sig[valid,k],
                            mu = hmm$mu[valid,k] + log(hmm$phi[m,k] * hmm$l[valid,k] * hmm$lambda[valid,k]),
                            log = TRUE
                        )
                } else {
                    l_y = 0
                }

                return(l_x + l_y)

            }
        ) %>% rowSums()

    })

    z = viterbi_compute(log(hmm$delta), logprob, hmm$logPi, N, M, nu, z)
        
    return(z)
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