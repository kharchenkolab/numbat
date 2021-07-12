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

Viterbi.dthmm.mv.inhom <- function (object, ...){
#     print('running multivariate inhomogenous HMM')
    x <- object$x
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    
    x2 <- object$x2
    dfunc2 <- HiddenMarkov:::makedensity(object$distn2)

    n <- length(x)
    m <- nrow(object$Pi[[1]])
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    
    nu[1, ] = log(object$delta)
        
    if (!is.na(x[1])) {
        nu[1, ] = nu[1, ] + log(dfunc(x=x[1], object$pm, HiddenMarkov:::getj(object$pn, 1)))
    }
    
    if (!is.na(x2[1])) {
        nu[1, ] = nu[1, ] + log(dfunc2(x=x2[1], object$pm2, HiddenMarkov:::getj(object$pn2, 1)))
    }
            
    logPi <- lapply(object$Pi, log)

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        nu[i, ] = apply(matrixnu + logPi[[i]], 2, max)
            
        if (!is.na(x[i])) {
            nu[i, ] = nu[i, ] + log(dfunc(x=x[i], object$pm, HiddenMarkov:::getj(object$pn, i)))
        }
        
        if (!is.na(x2[i])) {
            nu[i, ] = nu[i, ] + log(dfunc2(x=x2[i], object$pm2, HiddenMarkov:::getj(object$pn2, i)))
        }
    }
    
#     display(head(nu, 100))
#     display(dim(nu))

    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i+1]][, y[i+1]] + nu[i, ])
    return(y)
}


run_hmm_mv_inhom_norm = function(pAD, DP, p_s, logFC, sd, mu_neu, mu_del, mu_amp, p_0 = 1-1e-5, prior = NULL) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up",
               "5" = "loh_down", "6" = "amp_up", "7" = "amp_down")

    # intitial probabilities
    if (is.null(prior)) {
        w = 1/(1-p_0)
        prior = c((1 + w)/(length(states) + w), rep(1/(length(states) + w), length(states) - 1)) 
    }
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0, n_states) {
        matrix(
            c(
                p_0, rep((1-p_0)/6, 6),
                (1-p_0)/5, p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 4),
                (1-p_0)/5, p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 4),
                rep((1-p_0)/5, 3), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 3), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 5), p_0 * (1 - p_s), p_0 * p_s,
                rep((1-p_0)/5, 5), p_0 * p_s, p_0 * (1 - p_s)
            ),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0, n_states = length(states))}
    )
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(10, rep(c(10, 6), 3)), beta = c(10, rep(c(6, 10), 3))),
        pn = list(size = DP),
        discrete = TRUE
    )

    hmm$distn2 = 'norm'
    hmm$x2 = logFC
    hmm$pm2 = list(mean = c(mu_neu, rep(mu_del, 2), rep(mu_neu, 2), rep(mu_amp, 2)))
    hmm$pn2 = list(sd = sd)
    
    class(hmm) = 'dthmm.mv.inhom'

    return(states[as.character(HiddenMarkov::Viterbi(hmm))])
}

run_hmm_mv_inhom = function(pAD, DP, exp, p_s, sigma, mu_neu, mu_del, mu_gain, p_0 = 1-1e-5, prior = NULL) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up",
               "5" = "loh_down", "6" = "amp_up", "7" = "amp_down")

    # intitial probabilities
    if (is.null(prior)) {
        w = 1000
        prior = c((1 + w)/(length(states) + w), rep(1/(length(states) + w), length(states) - 1)) 
    }
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0, n_states) {
        matrix(
            c(
                p_0, rep((1-p_0)/6, 6),
                (1-p_0)/5, p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 4),
                (1-p_0)/5, p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 4),
                rep((1-p_0)/5, 3), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 3), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 5), p_0 * (1 - p_s), p_0 * p_s,
                rep((1-p_0)/5, 5), p_0 * p_s, p_0 * (1 - p_s)
            ),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0, n_states = length(states))}
    )
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(10, rep(c(10, 6), 3)), beta = c(10, rep(c(6, 10), 3))),
        pn = list(size = DP),
        discrete = TRUE
    )

    hmm$distn2 = 'norm'
    hmm$x2 = exp
    hmm$pm2 = list(mean = c(mu_neu, rep(mu_del, 2), rep(mu_neu, 2), rep(mu_gain, 2)), sd = rep(sigma, 7))
    
    class(hmm) = 'dthmm.mv.inhom'

    return(states[as.character(HiddenMarkov::Viterbi(hmm))])
}

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
    
    if (any(nu[n, ] == -Inf)) {
        stop("Problems With Underflow")
    }
#     display(head(mu, 100))
              
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i+1]][, y[i+1]] + nu[i, ])
        
    return(y)
}




run_hmm_mv_inhom_gpois = function(pAD, DP, p_s, Y_obs, lambda_ref, d_total, phi_neu, phi_del, phi_amp, alpha, beta, p_0 = 1-1e-5, gamma = 16, prior = NULL, exp_only = FALSE) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up",
               "5" = "loh_down", "6" = "amp_up", "7" = "amp_down")

    # intitial probabilities
    if (is.null(prior)) {
        w = 1/(1-p_0)
        prior = c((1 + w)/(length(states) + w), rep(1/(length(states) + w), length(states) - 1)) 
    }
    
    if (exp_only) {
        pAD = rep(NA, length(pAD))
        p_s = rep(0, length(p_s))
    }
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0, n_states) {
        matrix(
            c(
                p_0, rep((1-p_0)/6, 6),
                (1-p_0)/5, p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 4),
                (1-p_0)/5, p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 4),
                rep((1-p_0)/5, 3), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 3), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 5), p_0 * (1 - p_s), p_0 * p_s,
                rep((1-p_0)/5, 5), p_0 * p_s, p_0 * (1 - p_s)
            ),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0, n_states = length(states))}
    )
            
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(gamma * 0.5, rep(c(gamma * 0.6, gamma * 0.4), 3)), beta = c(gamma * 0.5, rep(c(gamma * 0.4, gamma * 0.6), 3))),
        pn = list(size = DP),
        discrete = TRUE
    )
    
    

    hmm$distn2 = 'gpois'
    hmm$x2 = Y_obs
    hmm$phi = c(phi_neu, rep(phi_del, 2), rep(phi_neu, 2), rep(phi_amp, 2))
    hmm$alpha = alpha
    hmm$beta = beta
    hmm$lambda_star = lambda_ref
    hmm$d = d_total
    
    class(hmm) = 'dthmm.mv.inhom.gpois'

    return(states[as.character(HiddenMarkov::Viterbi(hmm))])
}


run_hmm_mv_inhom_gpois2 = function(pAD, DP, p_s, Y_obs, lambda_ref, d_total, phi_neu, phi_del, phi_amp, phi_bamp, phi_bdel, alpha, beta, p_0 = 1-1e-5, gamma = 16, prior = NULL, exp_only = FALSE) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up",
               "5" = "loh_down", "6" = "amp_up", "7" = "amp_down", "8" = "bamp", "9" = "bdel")

    # intitial probabilities
    if (is.null(prior)) {
        w = 1/(1-p_0)
        prior = c((1 + w)/(length(states) + w), rep(1/(length(states) + w), length(states) - 1)) 
    }
    
    if (exp_only) {
        pAD = rep(NA, length(pAD))
        p_s = rep(0, length(p_s))
    }
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0, n_states) {
        matrix(
            c(
                p_0, rep((1-p_0)/8, 8),
                (1-p_0)/7, p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/7, 6),
                (1-p_0)/7, p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/7, 6),
                rep((1-p_0)/7, 3), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/7, 4),
                rep((1-p_0)/7, 3), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/7, 4),
                rep((1-p_0)/7, 5), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/7, 2),
                rep((1-p_0)/7, 5), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/7, 2),
                rep((1-p_0)/8, 7), p_0, (1-p_0)/8,
                rep((1-p_0)/8, 8), p_0
            ),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0, n_states = length(states))}
    )
            
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(gamma * 0.5, rep(c(gamma * 0.6, gamma * 0.4), 3), gamma * 0.5), beta = c(gamma * 0.5, rep(c(gamma * 0.4, gamma * 0.6), 3), gamma * 0.5, gamma * 0.5)),
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


Viterbi.dthmm.mv.inhom.poilog <- function (object, ...){

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
        
        nu[1, ] = nu[1, ] + sapply(
                object$pm2 * unlist(HiddenMarkov:::getj(object$pn2, 1)),
                function(mu) {
                    sads::dpoilog(x=x2[1], mu = log(mu), sig = object$sigma, log = TRUE)
                }
            )
        
#         mu[1, ] = c(
#             unlist(HiddenMarkov:::getj(object$pn2, 1)),
#             sapply(
#                 object$pm2 * unlist(HiddenMarkov:::getj(object$pn2, 1)),
#                 function(mu) {
#                     sads::dpoilog(x=x2[1], mu = log(mu), sig = object$overdisp, log = TRUE)
#                 }
#             ))
    }
            
    logPi <- lapply(object$Pi, log)

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        nu[i, ] = apply(matrixnu + logPi[[i]], 2, max)
            
        if (!is.na(x[i])) {
            nu[i, ] = nu[i, ] + dfunc(x=x[i], object$pm, HiddenMarkov:::getj(object$pn, i), log = TRUE)
        }
        
        if (!is.na(x2[i])) {
            nu[i, ] = nu[i, ] + sapply(
                object$pm2 * unlist(HiddenMarkov:::getj(object$pn2, i)),
                function(mu) {
                    sads::dpoilog(x=x2[i], mu = log(mu), sig = object$sigma, log = TRUE)
                }
            )
            
            
#             mu[i, ] = c(
#                 unlist(HiddenMarkov:::getj(object$pn2, i)),
#                 sapply(
#                 object$pm2 * unlist(HiddenMarkov:::getj(object$pn2, i)),
#                 function(mu) {
#                     sads::dpoilog(x=x2[i], mu = log(mu), sig = object$overdisp, log = TRUE)
#                 }
#             ))
        }
    }
    
#     display(head(mu))
#     display(head(nu))
              
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i+1]][, y[i+1]] + nu[i, ])
        
    return(y)
}

run_hmm_mv_inhom_poilog = function(pAD, DP, p_s, Y_obs, lambda_ref, d_total, phi_neu, phi_del, phi_amp, sigma, p_0 = 1-1e-5, prior = NULL, exp_only = FALSE) {
    
    # states
    states = c("1" = "neu", "2" = "del_up", "3" = "del_down", "4" = "loh_up",
               "5" = "loh_down", "6" = "amp_up", "7" = "amp_down")

    # intitial probabilities
    if (is.null(prior)) {
        w = 1/(1-p_0)
        prior = c((1 + w)/(length(states) + w), rep(1/(length(states) + w), length(states) - 1)) 
    }
    
    if (exp_only) {
        pAD = rep(NA, length(pAD))
        p_s = rep(0, length(p_s))
    }
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0, n_states) {
        matrix(
            c(
                p_0, rep((1-p_0)/6, 6),
                (1-p_0)/5, p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 4),
                (1-p_0)/5, p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 4),
                rep((1-p_0)/5, 3), p_0 * (1 - p_s), p_0 * p_s, rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 3), p_0 * p_s, p_0 * (1 - p_s), rep((1-p_0)/5, 2),
                rep((1-p_0)/5, 5), p_0 * (1 - p_s), p_0 * p_s,
                rep((1-p_0)/5, 5), p_0 * p_s, p_0 * (1 - p_s)
            ),
            ncol = n_states,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0, n_states = length(states))}
    )
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha = c(10, rep(c(10, 6), 3)), beta = c(10, rep(c(6, 10), 3))),
        pn = list(size = DP),
        discrete = TRUE
    )

    hmm$distn2 = 'poilog'
    hmm$x2 = Y_obs
    hmm$pm2 = c(phi_neu, rep(phi_del, 2), rep(phi_neu, 2), rep(phi_amp, 2))
    hmm$sigma = sigma
    hmm$pn2 = list(size = lambda_ref * d_total)
    
    class(hmm) = 'dthmm.mv.inhom.poilog'

    return(states[as.character(HiddenMarkov::Viterbi(hmm))])
}





