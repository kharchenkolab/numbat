
########## joint HMM ###########
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
run_hmm_mv_inhom = function(
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
                
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(
            alpha = alpha_states,
            beta = beta_states
        ),
        pn = DP,
        discrete = TRUE
    )

    hmm$x2 = Y_obs
    hmm$phi = phi_states
    hmm$lambda_star = lambda_ref
    hmm$d = d_total

    if (length(mu) == 1 & length(sig) == 1) {
        mu = rep(mu, length(Y_obs))
        sig = rep(sig, length(Y_obs))
    }

    hmm$distn2 = 'poilog'
    hmm$mu = mu
    hmm$sig = sig
    
    class(hmm) = 'dthmm.mv.inhom.lnpois'

    MPC = viterbi_joint(hmm)
        
    return(states[as.character(MPC$y)])
}


run_hmm_mv_inhom2 = function(
    pAD_1, pAD_2, DP_1, DP_2, Y_obs_1, Y_obs_2, d_total_1, d_total_2, lambda_ref_1, lambda_ref_2, p_s, theta_min = 0.08, theta_neu = 0,
    bal_cnv = TRUE, phi_del = 2^(-0.25), phi_amp = 2^(0.25), phi_bamp = phi_amp, phi_bdel = phi_del, 
    alpha = 1, beta = 1, 
    mu_1 = 0, sig_1 = 1, mu_2 = 0, sig_2 = 1,
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

    N = length(Y_obs_1)

    if (length(mu_1) == 1 & length(sig_1) == 1) {
        mu_1 = rep(mu_1, N)
        sig_1 = rep(sig_1, N)
    }

    if (length(mu_2) == 1 & length(sig_2) == 1) {
        mu_2 = rep(mu_2, N)
        sig_2 = rep(sig_2, N)
    }

    if (length(d_total_1) == 1 & length(d_total_2) == 1) {
        d_total_1 = rep(d_total_1, N)
        d_total_2 = rep(d_total_2, N)
    }
                
    hmm = list(
        x = matrix(c(pAD_1, pAD_2), ncol = 2), 
        d = matrix(c(DP_1, DP_2), ncol = 2),
        y = matrix(c(Y_obs_1, Y_obs_2), ncol = 2),
        l = matrix(c(d_total_1, d_total_2), ncol = 2),
        lambda = matrix(c(lambda_ref_1, lambda_ref_2), ncol = 2),
        mu = matrix(c(mu_1, mu_2), ncol = 2),
        sig = matrix(c(sig_1, sig_2), ncol = 2),
        logPi = log(As) * 2, 
        phi = matrix(rep(phi_states, 2), ncol = 2),
        delta = matrix(rep(prior, 2), ncol = 2), 
        alpha = matrix(rep(alpha_states, 2), ncol = 2),
        beta = matrix(rep(beta_states, 2), ncol = 2),
        states = matrix(rep(states, 2), ncol = 2)
    )

    MPC = viterbi_joint2(hmm)
        
    return(tibble(state = states[as.character(MPC$z)], max_LL = MPC$LL))
}

############ allele HMM ############

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

#' @keywords internal  
makedensity <- function(distn){
    ## https://github.com/cran/HiddenMarkov/blob/master/R/makedensity.R
    dname <- paste("d", distn[1], sep="")
    x <- paste("function(x, pm, pn=NULL, log=FALSE)
         do.call(\"", dname, "\", c(list(x=x), pm, pn,", sep="")
    if (distn[1]=="glm") x <- paste(x, " list(family=\"", distn[2],
         "\", link=\"", distn[3], "\"),", sep="")
    eval(parse(text=paste(x, " list(log=log)))", sep="")))
}


#' @keywords internal
getj <- function(x, j){
    ## https://github.com/cran/HiddenMarkov/blob/master/R/getj.R
    #   get the jth "row" from a list
    if (is.null(x)) return(NULL)
    n <- length(x)
    for (i in 1:n)
        x[[i]] <- x[[i]][j]
    return(x)
}

############ time inhomogenous univariate HMM ############

#' Viterbi algorithm for allele HMM
#' @keywords internal
viterbi_allele <- function(hmm) {

    N <- length(hmm$x)
    M <- nrow(hmm$Pi[[1]])
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)

    logPi <- lapply(hmm$Pi, log)

    logprob = sapply(1:M, function(m) {

        l_x = dbbinom(x = hmm$x, size = hmm$d, alpha = hmm$alpha[m], beta = hmm$beta[m], log = TRUE)

        l_x[is.na(l_x)] = 0

        return(l_x)

    })

    nu[1, ] <- log(hmm$delta) + logprob[1,]

    for (i in 2:N) {
        matrixnu <- matrix(nu[i - 1, ], nrow = M, ncol = M)
        nu[i, ] = apply(matrixnu + logPi[[i]], 2, max)
        nu[i, ] = nu[i, ] + logprob[i,]
    }

    z[N] <- which.max(nu[N, ])

    for (i in seq(N - 1, 1, -1)) z[i] <- which.max(logPi[[i + 1]][, z[i + 1]] + nu[i, ])

    LL = max(nu[N, ])
        
    return(z)
}

#' Forward-backward algorithm for allele HMM
#' @keywords internal
forward_back_allele = function (obj, ...) {

    # case of one-data point
    if (length(obj$x) == 1) {
        return(NA)
    }
    
    x <- obj$x
    p_x <- makedensity(obj$distn)
    
    m <- nrow(obj$Pi[[1]])
    n <- length(x)
    
    logprob = sapply(1:m, function(k) {
        
        l_x = p_x(x = x,  getj(obj$pm, k), list(size = obj$pn), log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })
        
    logphi <- log(as.double(obj$delta))
    
    logPi <- lapply(obj$Pi, log)
        
    marginals = forward_backward_compute(obj, logphi, logprob, logPi, n, m)

    colnames(marginals) = obj$states

    return(marginals)
}


# Only compute total log likelihood
#' @keywords internal
likelihood_allele = function (obj, ...) {
        
    x <- obj$x
    p_x <- makedensity(obj$distn)
    
    m <- nrow(obj$Pi[[1]])
    n <- length(x)
    
    logprob = sapply(1:m, function(k) {
        
        l_x = p_x(x = x,  getj(obj$pm, k), list(size = obj$pn), log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })

    logphi <- log(as.double(obj$delta))
    
    logPi <- lapply(obj$Pi, log)
        
    LL <- likelihood_allele_compute(obj, logphi, logprob, logPi, n, m)

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
run_hmm_inhom = function(pAD, DP, p_s, t = 1e-5, theta_min = 0.08, gamma = 20, prior = NULL) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('More than one gamma parameter')
    }
    
    # states
    states = c("neu", "theta_1_up", "theta_1_down", "theta_2_up", "theta_2_down")

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
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, t)}
    )

    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/5, 5)
    }

    theta_1 = theta_min
    theta_2 = 0.4
            
    hmm = list(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        alpha = gamma * c(0.5, 0.5 + theta_1, 0.5 - theta_1, 0.5 + theta_2, 0.5 - theta_2),
        beta = gamma * c(0.5, 0.5 - theta_1, 0.5 + theta_1, 0.5 - theta_2, 0.5 + theta_2),
        d = DP
    )
    
    solution = states[viterbi_allele(hmm)]
    
    return(solution)
}


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
    calc_trans_mat = function(p_s) {
        matrix(c(1 - p_s, p_s, p_s, 1 - p_s), ncol = 2, byrow = TRUE)
    }
    As = lapply(p_s, function(p_s) {
        calc_trans_mat(p_s)
    })
    prior = c(0.5, 0.5)
    alpha_up = (0.5 + theta) * gamma
    beta_up = (0.5 - theta) * gamma
    alpha_down = beta_up
    beta_down = alpha_up
    
    hmm = HiddenMarkov::dthmm(x = pAD, Pi = As, delta = prior, 
        distn = "bbinom", pm = list(alpha = c(alpha_up, alpha_down), 
            beta = c(beta_up, beta_down)), pn = DP, 
        discrete = TRUE)

    hmm$states = states

    class(hmm) = "dthmm.inhom"

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

############ Joint HMM ############

#' Viterbi algorithm for joint HMM (PLN expression likelihood)
#' @keywords internal
viterbi_joint <- function (object, ...){

    x <- object$x
    dfunc <- makedensity(object$distn)
    
    x2 <- object$x2
    dfunc2 <- makedensity(object$distn2)

    n <- length(x)
    m <- nrow(object$Pi[,,1])
    nu <- matrix(NA, nrow = n, ncol = m)
    mu <- matrix(NA, nrow = n, ncol = m + 1)
    y <- rep(NA, n)
    
    nu[1, ] = log(object$delta)
    
        
    if (!is.na(x[1])) {
        nu[1, ] = nu[1, ] + dfunc(x=x[1], object$pm, object$pn[1], log = TRUE)
    }
    
    if (!is.na(x2[1])) {

        nu[1, ] = nu[1, ] + dfunc2(
            x = rep(x2[1], m),
            list('sig' = rep(object$sig[1], m)),
            list('mu' = object$mu[1] + log(object$phi * object$d * object$lambda_star[1])),
            log = TRUE
        )
        
    }
            
    logPi <- log(object$Pi)

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        nu[i, ] = apply(matrixnu + logPi[,,i], 2, max)
            
        if (!is.na(x[i])) {
            nu[i, ] = nu[i, ] + dfunc(x=x[i], object$pm, object$pn[i], log = TRUE)
        }
        
        if (!is.na(x2[i])) {
            nu[i, ] = nu[i, ] + dfunc2(
                x = rep(x2[i], m),
                list('sig' = rep(object$sig[i], m)),
                list('mu' = object$mu[i] + log(object$phi * object$d * object$lambda_star[i])),
                log = TRUE
            )
        }
    }

    if (any(is.na(nu))) {
        # fwrite(nu, '~/debug.txt')
        stop("NA values in viterbi")
    }
    
    if (all(nu[n, ] == -Inf)) {
        # fwrite(nu, '~/debug.txt')
        stop("Problems With Underflow")
    }

    # fwrite(nu, '~/debug.txt')
              
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[,,i+1][, y[i+1]] + nu[i, ])

    LL = max(nu[n, ])

    # print(LL)
        
    return(list(y = y, LL = LL))
}

#' Viterbi for clonal LOH detection
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


#' Calculate the transition matrix for joint HMM
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
# cn/phase are sclars, only p_s is vectorized
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

get_joint_hmm = function(
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
    logPi_cn = calc_trans_mat2(t, p_s, w, states_cn, states_phase, part = 'cn')
    logPi_phase = calc_trans_mat2(t, p_s, w, states_cn, states_phase, part = 'phase')
    logPi = logPi_cn + logPi_phase
    # logPi = calc_trans_mat(t, p_s, w, states_cn, states_phase) %>% log

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
    logPi = logPi[states_index, states_index,]
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
        logPi = logPi, 
        logPi_cn = logPi_cn, 
        logPi_phase = logPi_phase, 
        phi = matrix(phi_states),
        delta = matrix(prior), 
        alpha = matrix(alpha_states),
        beta = matrix(beta_states),
        states = matrix(states),
        p_s = p_s
    )

    return(hmm)
}


#' Viterbi algorithm for joint HMM (PLN expression likelihood)
#' @keywords internal
viterbi_joint2 <- function (hmm, ...){

    N <- nrow(hmm$x)
    M <- nrow(hmm$logPi[,,1])
    K <- ncol(hmm$x)
    nu <- matrix(NA, nrow = N, ncol = M)
    z <- rep(NA, N)
    
    nu[1, ] = rowSums(log(hmm$delta))

    logprob = sapply(1:M, function(m) {

        sapply(
            1:K,
            function(k) {

                l_x = dbbinom(x = hmm$x[,k], size = hmm$d[,k], alpha = hmm$alpha[m,k], beta = hmm$beta[m,k], log = TRUE)

                l_x[is.na(l_x)] = 0

                l_y = rep(0,N)
                valid = !is.na(hmm$y[,k])

                l_y[valid] = dpoilog(
                        x = hmm$y[valid,k],
                        sig = hmm$sig[valid,k],
                        mu = hmm$mu[valid,k] + log(hmm$phi[m,k] * hmm$l[valid,k] * hmm$lambda[valid,k]),
                        log = TRUE
                    )

                return(l_x + l_y)

            }
        ) %>% rowSums()

    })

    for (i in 1:N) {

        if (i > 1) {
            matrixnu <- matrix(nu[i - 1, ], nrow = M, ncol = M)
            nu[i, ] = apply(matrixnu + hmm$logPi[,,i], 2, max)
        }

        nu[i, ] = nu[i, ] + logprob[i, ]
    }

    if (any(is.na(nu))) {
        stop("NA values in viterbi")
    }
    
    if (all(nu[n, ] == -Inf)) {
        stop("Problems With Underflow")
    }
              
    z[N] <- which.max(nu[N, ])

    for (i in seq(N - 1, 1, -1)) z[i] <- which.max(hmm$logPi[,,i+1][, z[i+1]] + nu[i, ])

    LL = max(nu[N, ])
        
    return(list(z = hmm$states[z,], LL = LL))
}

#' Forward-backward algorithm for joint HMM
#' @keywords internal
likelihood_joint = function(hmm) {
    
    M <- nrow(hmm$logPi[,,1])
    N <- nrow(hmm$x)
    K = ncol(hmm$x)

    logPi = lapply(1:N, function(i) hmm$logPi[,,i])

    # case of one-data point
    if (N == 1) {
        return(NA)
    }
    
    logprob = sapply(1:M, function(m) {

        sapply(
            1:K,
            function(k) {

                l_x = dbbinom(x = hmm$x[,k], size = hmm$d[,k], alpha = hmm$alpha[m,k], beta = hmm$beta[m,k], log = TRUE)

                l_x[is.na(l_x)] = 0

                l_y = rep(0,N)
                valid = !is.na(hmm$y[,k])

                l_y[valid] = dpoilog(
                        x = hmm$y[valid,k],
                        sig = hmm$sig[valid,k],
                        mu = hmm$mu[valid,k] + log(hmm$phi[m,k] * hmm$l[valid,k] * hmm$lambda[valid,k]),
                        log = TRUE
                    )

                return(l_x + l_y)

            }
        ) %>% rowSums()

    })
        
    logphi <- rowSums(log(hmm$delta))
            
    LL = likelihood_allele_compute(hmm, logphi, logprob, logPi, N, M)

    return(LL)
}

# bundle two HMMs using different transition kernels
merge_hmm = function(hmm1, hmm2, kernel = 'indep') {
    
    N = dim(hmm1$logPi)[3]
    M = dim(hmm1$logPi)[2]
    
    states = cross(hmm1$states, hmm2$states)
    
    states_cn = states %>% str_remove('_up|_down') %>% matrix(ncol = 2)
    states_ph = states %>% str_extract('up|down') %>% matrix(ncol = 2)

    if (kernel == 'equal') {

        keep = which(states[,1] == states[,2])
        logPi = hmm1$logPi + hmm2$logPi_cn
        logPi_ph = hmm1$logPi_ph

    } else if (kernel == 'indep') {

        keep = 1:nrow(states)

        logPi = sapply(
            1:N,
            function(n) {
                cross_trans_indep(hmm1$logPi[,,n], hmm2$logPi[,,n])
            }
        ) %>%
        array(dim = c(M^2, M^2, N))

    } else if (kernel == 'phase') {

        keep = 1:nrow(states)
        
        logPi = calc_trans_mat3(t, hmm1$p_s, w, states_cn, states_ph, kernel = 'phase')

    }
    
    list(
        states = states[keep,],
        x = matrix(c(hmm1$x, hmm2$x), ncol = 2),
        d = matrix(c(hmm1$d, hmm2$d), ncol = 2),
        y = matrix(c(hmm1$y, hmm2$y), ncol = 2),
        l = matrix(c(hmm1$l, hmm2$l), ncol = 2),
        lambda = matrix(c(hmm1$lambda, hmm2$lambda), ncol = 2),
        mu = matrix(c(hmm1$mu, hmm2$mu), ncol = 2),
        sig = matrix(c(hmm1$sig, hmm2$sig), ncol = 2),
        phi = cross(hmm1$phi, hmm2$phi)[keep,],
        alpha = cross(hmm1$alpha, hmm2$alpha)[keep,],
        beta = cross(hmm1$beta, hmm2$beta)[keep,],
        delta = cross(hmm1$delta, hmm2$delta)[keep,],
        logPi = logPi,
        logPi_ph = logPi_ph
    )
}

cross = function(x,y) {
    m = as.matrix(expand.grid(x, y))
    colnames(m) = c('z1', 'z2')
    return(m)
}

cross_trans_indep = function(logA1, logA2, opt = '+') {
    
    m = nrow(logA1)
    array(aperm(outer(logA1, logA2, opt), c(1,3,2,4)), dim = c(m^2, m^2))

}

renorm = function(logA) {
    t(apply(logA, 1, function(x){x - logSumExp(x)}))
}

compare_hmm = function(hmm1, hmm2) {
    hmm_merged = merge_hmm(hmm1, hmm2, kernel = 'phase')
    L1 = viterbi_joint2(hmm_merged)$LL
    hmm_merged = merge_hmm(hmm1, hmm2, kernel = 'equal')
    L0 = viterbi_joint2(hmm_merged)$LL
    return(L1 - L0)
}

calc_trans_mat2 = function(t, p_s, w, states_cn, states_phase, part = 'cn') {
    
    if (part == 'cn') {
        get_trans_probs = get_trans_probs_cn
    } else {
        get_trans_probs = get_trans_probs_phase
    }

    sapply(1:length(states_cn), function(from) {
        sapply(1:length(states_cn), function(to) {
            get_trans_probs(t, p_s, w, states_cn[from], states_phase[from], states_cn[to], states_phase[to])
        }) %>% t
    }) %>% t %>%
    array(dim = c(length(states_cn), length(states_cn), length(p_s)))

}

get_trans_probs_cn = function(t, p_s, w, cn_from, phase_from, cn_to, phase_to) {

    if (cn_from == cn_to) {
        logp = rep(log(1-t), length(p_s))
    } else {
        logp = log(t) + log(w[[cn_to]]/sum(w[names(w)!=cn_from]))
        if (!is.na(phase_to)) {
            logp = logp - log(2)
        }
        logp = rep(logp, length(p_s))
    }
    
    return(logp)
}

get_trans_probs_phase = function(t, p_s, w, cn_from, phase_from, cn_to, phase_to) {

    # if (cn_from == cn_to) {
    #     if (is.na(phase_from) & is.na(phase_to)) {
    #         logp = rep(0, length(p_s))
    #     } else if (phase_from == phase_to) {
    #         logp = log(1-p_s)
    #     } else {
    #         logp = log(p_s)
    #     }
    # } else {
    #     logp = rep(0, length(p_s))
    # }

    if (is.na(phase_from) | is.na(phase_to)) {
        logp = rep(0, length(p_s))
    } else if (phase_from == phase_to) {
        logp = log(1-p_s)
    } else {
        logp = log(p_s)
    }
    
    return(logp)
}

calc_trans_mat3 = function(t, p_s, w, states_cn, states_ph, kernel = 'phase') {

    N = length(p_s)
    M = nrow(states_cn)

    sapply(1:M, function(from) {
        sapply(1:M, function(to) {
            get_trans_probs_cn3(t, w, states_cn[from,1], states_cn[to,1], N) +
            get_trans_probs_cn3(t, w, states_cn[from,2], states_cn[to,2], N) +
            get_trans_probs_ph3(t, p_s, states_ph[from,], states_ph[to,], N, kernel = kernel)
        }) %>% t
    }) %>% t %>%
    array(dim = c(M, M, N))

}

get_trans_probs_cn3 = function(t, w, cn_from, cn_to, N) {

    if (cn_from == cn_to) {
        logp = rep(log(1-t), N)
    } else {
        logp = log(t) + log(w[[cn_to]]/sum(w[names(w)!=cn_from]))
        if (!cn_to %in% c('neu', 'bamp', 'bdel')) {
            logp = logp - log(2)
        }
        logp = rep(logp, N)
    }
    
    return(logp)
}

get_trans_probs_ph3 = function(t, p_s, phase_from, phase_to, N, kernel = 'phase') {

    switch = phase_from != phase_to

    if (is.na(switch[1]) & is.na(switch[2])) {
        logp = rep(0, N)
    } else if (is.na(switch[1])) {
        if (switch[2]) {
            logp = log(p_s)
        } else {
            logp = log(1-p_s)
        }
    } else if (is.na(switch[2])) {
        if (switch[1]) {
            logp = log(p_s)
        } else {
            logp = log(1-p_s)
        }
    } else if ((switch[1] & !switch[2]) | (!switch[1] & switch[2])) {
        if (kernel == 'phase') {
            logp = rep(-Inf, N)
        } else {
            logp = log(p_s)
        }
    } else if ((!switch[1]) & (!switch[2])) {
        logp = log(1-p_s)
    } else if (switch[1] & switch[2]) {
        logp = log(p_s)
    }

    return(logp)
}


