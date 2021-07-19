require(extraDistr)

############ time homogenous univariate HMM ############

run_hmm = function(pAD, DP, p_0 = 1-1e-5, p_s = 0.1) {
    
    # states
    states = c("theta_up", "neu", "theta_down")
    
    # transition matrix
    A <- matrix(
        c(p_0 * (1-p_s), 1 - p_0, p_0 * p_s, 
         (1-p_0)/2, p_0, (1-p_0)/2,
         p_0 * p_s, 1-p_0, p_0 * (1 - p_s)),
        ncol = length(states),
        byrow = TRUE
    )

    # intitial probabilities
    prior = rep(1/length(states), length(states))
    
    hmm = HiddenMarkov::dthmm(
            x = pAD, 
            Pi = A, 
            delta = prior, 
            distn = "bbinom",
            pm = list(alpha=c(10,10,6), beta=c(6,10,10)),
            pn = list(size = DP),
            discrete=TRUE)
    
    solution = states[HiddenMarkov::Viterbi(hmm)]
    
    return(solution)
}

############ time inhomogenous univariate HMM ############

# phase switch probablity as a function of genomic distance
switch_prob = function(x, pad = 0.01) {
    1-pmax(pmin(2.8 - 0.38 * log10(x), 1 - pad), 0.5 + pad)
}

Viterbi.dthmm.inhom <- function (object, ...){
#     print('Solving univariate nonhomogenous markov chain')
    x <- object$x
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    n <- length(x)
    m <- nrow(object$Pi[[1]])
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    nu[1, ] <- log(object$delta) + dfunc(x=x[1], object$pm,
                                    HiddenMarkov:::getj(object$pn, 1), log=TRUE)
    logPi <- lapply(object$Pi, log)
    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi[[i]], 2, max) +
            dfunc(x=x[i], object$pm, HiddenMarkov:::getj(object$pn, i), log=TRUE)
    }
    
#     if (any(nu[n, ] == -Inf)) 
#         stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    # double check this index of logPi
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i + 1]][, y[i + 1]] + nu[i, ])
        
    max.loglik = max(nu[n, ])
    
    return(list('y' = y, 'max.loglik' = max.loglik))
}

run_hmm_inhom = function(pAD, DP, p_s, p_0 = 1-1e-5, prior = NULL) {
    
    # states
    states = c("theta_up", "neu", "theta_down")
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0, n_states) {
        matrix(
            c(p_0 * (1-p_s), 1 - p_0, p_0 * p_s, 
             (1-p_0)/2, p_0, (1-p_0)/2,
             p_0 * p_s, 1-p_0, p_0 * (1 - p_s)),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0, n_states = length(states))}
    )
    
    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/length(states), length(states))
    }
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha=c(10,10,6), beta=c(6,10,10)),
        pn = list(size = DP),
        discrete = TRUE)

    class(hmm) = 'dthmm.inhom'
        
    solution = states[HiddenMarkov::Viterbi(hmm)$y]
    
    return(solution)
}


############ time homogenous multivariate HMM ############

Viterbi.dthmm.mv <- function (object, ...){
#     print('running multivariate HMM')
    x <- object$x
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    
    x2 <- object$x2
    dfunc2 <- HiddenMarkov:::makedensity(object$distn2)
    
    n <- length(x)
    m <- nrow(object$Pi)
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    
    nu[1, ] = log(object$delta)
    
    if (!is.na(x[1])) {
        nu[1, ] = nu[1, ] + dfunc(x=x[1], object$pm, HiddenMarkov:::getj(object$pn, 1), log=TRUE)
    }
    
    if (!is.na(x2[1])) {
        nu[1, ] = nu[1, ] + dfunc2(x=x2[1], object$pm2, HiddenMarkov:::getj(object$pn2, 1), log=TRUE)
    }
        
    logPi <- log(object$Pi)

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        nu[i, ] = apply(matrixnu + logPi, 2, max)
            
        if (!is.na(x[i])) {
            nu[i, ] = nu[i, ] + dfunc(x=x[i], object$pm, HiddenMarkov:::getj(object$pn, i), log=TRUE)
        }
        
        if (!is.na(x2[i])) {
            nu[i, ] = nu[i, ] + dfunc2(x=x2[i], object$pm2, HiddenMarkov:::getj(object$pn2, i), log=TRUE)
        }
    }
    
    if (any(nu[n, ] == -Inf)) 
        stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])
    return(y)
}

run_hmm_mv = function(pAD, DP, exp, sigma, mu_neu, mu_del, mu_gain, p_0 = 1-1e-5, p_s = 0.1) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up",
               "5" = "loh_down", "6" = "amp_up", "7" = "amp_down")

    # intitial probabilities
    prior = rep(1/length(states), length(states))

    # transition matrix
    A <- matrix(
        c(
            p_0, rep((1-p_0)/6, 6),
            (1-p_0)/5, p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 4),
            (1-p_0)/5, p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 4),
            rep((1-p_0)/5, 3), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 2),
            rep((1-p_0)/5, 3), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 2),
            rep((1-p_0)/5, 5), p_0 * (1 - p_s), p_0 * p_s,
            rep((1-p_0)/5, 5), p_0 * p_s, p_0 * (1 - p_s)
        ),
        ncol = length(states),
        byrow = TRUE
    )
    
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = A, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(10, rep(c(10, 6), 3)), beta = c(10, rep(c(6, 10), 3))),
        pn = list(size = DP),
        discrete = TRUE
    )

    hmm$distn2 = 'norm'
    hmm$x2 = exp
    hmm$pm2 = list(mean = c(mu_neu, rep(mu_del, 2), rep(mu_neu, 2), rep(mu_gain, 2)), sd = rep(sigma, 7))
    
    class(hmm) = 'dthmm.mv'

    return(states[as.character(HiddenMarkov::Viterbi(hmm))])
}

############ time inhomogenous multivariate HMM ############
Viterbi.dthmm.mv.inhom.gpois <- function (object, ...){

    x <- object$x
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    
    x2 <- object$x2
    dfunc2 <- HiddenMarkov:::makedensity(object$distn2)

    n <- length(x)
    m <- nrow(object$Pi[[1]])
    nu <- matrix(NA, nrow = n, ncol = m)
    mu <- matrix(NA, nrow = n, ncol = m + 1)
    y <- rep(NA, n)
    
    
    nu[1, ] = log(object$delta)
    
        
    if (!is.na(x[1])) {
        nu[1, ] = nu[1, ] + dfunc(x=x[1], object$pm, HiddenMarkov:::getj(object$pn, 1), log = TRUE)
    }
    
    if (!is.na(x2[1])) {

        nu[1, ] = nu[1, ] + dfunc2(
            x = x2[1],
            list('shape' = object$alpha),
            list('rate' = object$beta/(object$phi * object$d * object$lambda_star[1])),
            log = TRUE
        )
        
    }
            
    logPi <- lapply(object$Pi, log)

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        nu[i, ] = apply(matrixnu + logPi[[i]], 2, max)
            
        if (!is.na(x[i])) {
            nu[i, ] = nu[i, ] + dfunc(x=x[i], object$pm, HiddenMarkov:::getj(object$pn, i), log = TRUE)
        }
        
        if (!is.na(x2[i])) {
            nu[i, ] = nu[i, ] + dfunc2(
                x = x2[i],
                list('shape' = object$alpha),
                list('rate' = object$beta/(object$phi * object$d * object$lambda_star[i])),
                log = TRUE
            )
            
        }
    }
    
#     if (any(nu[n, ] == -Inf)) {
#         stop("Problems With Underflow")
#     }
#     display(head(mu, 100))
              
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i+1]][, y[i+1]] + nu[i, ])
        
    return(y)
}

get_trans_probs = function(t, w) {
    
    a = list()
        
    for (from in names(w)) {
        for (to in names(w)) {
            if (from == to) {
                a[[from]][[to]] = 1-t
            } else {
                a[[from]][[to]] = t * w[[to]]/sum(w[names(w)!=from])
            }
        }
    }
    
    return(a)
}

run_hmm_mv_inhom_gpois = function(pAD, DP, p_s, Y_obs, lambda_ref, d_total, bal_cnv = TRUE, phi_neu = 1, phi_del = 2^(-0.25), phi_amp = 2^(0.25), phi_bamp = 2^(0.5), phi_bdel = 2^(-0.5), alpha = 1, beta = 1, t = 1e-5, gamma = 18, prior = NULL, exp_only = FALSE, allele_only = FALSE) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up", "5" = "loh_down", 
               "6" = "amp_up", "7" = "amp_down", "8" = "bamp", "9" = "bdel")
    
    # relative abundance of states
    w = c('neu' = 1, 'del' = 1, 'loh' = 1, 'amp' = 1, 'bamp' = 1e-5, 'bdel' = 1e-10)
    
    if (!bal_cnv) {
        w[c('bamp', 'bdel')] = 0
    }

    # intitial probabilities
    if (is.null(prior)) {
        # encourage CNV from telomeres
        a_0 = get_trans_probs(t = t * 100, w)[['neu']]
        prior = c(a_0[['neu']], 
            rep(a_0[['del']]/2, 2),
            rep(a_0[['loh']]/2, 2),
            rep(a_0[['amp']]/2, 2), 
            a_0[['bamp']],
            a_0[['bdel']])
    }
    
    if (exp_only) {
        pAD = rep(NA, length(pAD))
        p_s = rep(0, length(p_s))
    }
    
    if (allele_only) {
        Y_obs = rep(NA, length(Y_obs))
    }
    
    # transition matrices
    calc_trans_mat = function(p_s, t, n_states) {
        a = get_trans_probs(t, w)
        matrix(
            c(
                1-t, rep(a[['neu']][['del']]/2, 2), rep(a[['neu']][['loh']]/2, 2), rep(a[['neu']][['amp']]/2, 2), a[['neu']][['bamp']], a[['neu']][['bdel']],
                a[['del']][['neu']], (1-t)*(1-p_s), (1-t)*p_s, rep(a[['del']][['loh']]/2, 2), rep(a[['del']][['amp']]/2, 2), a[['del']][['bamp']], a[['del']][['bdel']],
                a[['del']][['neu']], (1-t)*p_s, (1-t)*(1-p_s), rep(a[['del']][['loh']]/2, 2), rep(a[['del']][['amp']]/2, 2), a[['del']][['bamp']], a[['del']][['bdel']],
                a[['loh']][['neu']], rep(a[['loh']][['del']]/2, 2), (1-t)*(1-p_s), (1-t)*p_s, rep(a[['loh']][['amp']]/2, 2), a[['loh']][['bamp']], a[['loh']][['bdel']],
                a[['loh']][['neu']], rep(a[['loh']][['del']]/2, 2), (1-t)*p_s, (1-t)*(1-p_s), rep(a[['loh']][['amp']]/2, 2), a[['loh']][['bamp']], a[['loh']][['bdel']],
                a[['amp']][['neu']], rep(a[['amp']][['del']]/2, 2), rep(a[['amp']][['loh']]/2, 2), (1-t)*(1-p_s), (1-t)*p_s, a[['amp']][['bamp']], a[['amp']][['bdel']],
                a[['amp']][['neu']], rep(a[['amp']][['del']]/2, 2), rep(a[['amp']][['loh']]/2, 2), (1-t)*p_s, (1-t)*(1-p_s), a[['amp']][['bamp']], a[['amp']][['bdel']],
                a[['bamp']][['neu']], rep(a[['bamp']][['del']]/2, 2), rep(a[['bamp']][['loh']]/2, 2), rep(a[['bamp']][['amp']]/2, 2), 1-t, a[['bamp']][['bdel']],
                a[['bdel']][['neu']], rep(a[['bdel']][['del']]/2, 2), rep(a[['bdel']][['loh']]/2, 2), rep(a[['bdel']][['amp']]/2, 2), a[['bdel']][['bamp']], 1-t
            ),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, t, n_states = length(states))}
    )
    
    theta_u = 0.58
    theta_d = 0.42
            
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(gamma * 0.5, rep(c(gamma * theta_u, gamma * theta_d), 3), gamma * 0.5), beta = c(gamma * 0.5, rep(c(gamma * theta_d, gamma * theta_u), 3), gamma * 0.5, gamma * 0.5)),
        pn = list(size = DP),
        discrete = TRUE
    )
    
    hmm$distn2 = 'gpois'
    hmm$x2 = Y_obs
    hmm$phi = c(phi_neu, rep(phi_del, 2), rep(phi_neu, 2), rep(phi_amp, 2), phi_bamp, phi_bdel)
    hmm$alpha = alpha
    hmm$beta = beta
    hmm$lambda_star = lambda_ref
    hmm$d = d_total
    
    class(hmm) = 'dthmm.mv.inhom.gpois'

    return(states[as.character(HiddenMarkov::Viterbi(hmm))])
}






