
########## joint HMM ###########
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
        pn = list(size = DP),
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

############ allele HMM ############

#' Beta-binomial distribution density function
#' A distribution is beta-binomial if p, the probability of success, 
#' in a binomial distribution has a beta distribution with shape 
#' parameters α > 0 and β > 0
#' For more details, see extraDistr::dbbinom
#'
#' @param x vector of quantiles
#' @param size number of trials (zero or more)
#' @param alpha numeric (default=1)
#' @param beta numeric (default=1)
#' @param log boolean (default=FALSE)
#' @return density values returned as numeric vector
#' @examples
#' xx <- 1:1000
#' dbbinom(xx, 1000, 5, 13)
#'
#' @export
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
viterbi_allele <- function (obj, ...){
#     print('Solving univariate nonhomogenous markov chain')
    x <- obj$x
    dfunc <- makedensity(obj$distn)
    n <- length(x)
    m <- nrow(obj$Pi[[1]])
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)

    nu[1, ] <- log(obj$delta) 

    if (!is.na(x[1])) {
        nu[1, ] <- nu[1, ] + dfunc(x=x[1], obj$pm, getj(obj$pn, 1), log=TRUE)
    }

    if (is.null(obj$logPi)) {
        logPi <- lapply(obj$Pi, log)
    } else {
        logPi = obj$logPi
    }

    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi[[i]], 2, max)
        if (!is.na(x[i])) {
            nu[i, ] <- nu[i, ] + dfunc(x=x[i], obj$pm, getj(obj$pn, i), log=TRUE)
        }
    }
#     if (any(nu[n, ] == -Inf)) 
#         stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    # double check this index of logPi
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[[i + 1]][, y[i + 1]] + nu[i, ])

    LL = max(nu[n, ])
        
    return(y)
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

    hmm = get_full_allele_hmm(pAD, DP, p_s, t = t, theta_min = theta_min, gamma = gamma, prior = prior)
        
    solution = hmm$states[viterbi_allele(hmm)]
    
    return(solution)
}

#' Allele HMM with two theta states
#' @keywords internal
get_full_allele_hmm = function(pAD, DP, p_s, t = 1e-5, theta_min = 0.08, gamma = 20, prior = NULL) {

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
            
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(
            alpha = gamma * c(0.5, 0.5 + theta_1, 0.5 - theta_1, 0.5 + theta_2, 0.5 - theta_2),
            beta = gamma * c(0.5, 0.5 - theta_1, 0.5 + theta_1, 0.5 - theta_2, 0.5 + theta_2)
        ),
        pn = list(size = DP),
        discrete = TRUE
    )

    hmm$states = states

    class(hmm) = 'dthmm.inhom'

    return(hmm)
}

#' Allele HMM with one theta state
#' @keywords internal
get_single_allele_hmm = function(pAD, DP, p_s, t = 1e-5, theta_min = 0.08, gamma = 20, theta_neu = 0, prior = NULL) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('More than one gamma parameter')
    }
    
    # states
    states = c("neu", "theta_up", "theta_down")
    
    # transition matrices
    calc_trans_mat = function(p_s, t) {
        matrix(
            c(
                (1-t), t/2, t/2,
                t, (1-t)*(1-p_s), (1-t)*p_s, 
                t, (1-t)*p_s, (1-t)*(1-p_s)
            ),
            ncol = 3,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, t)}
    )
    
    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/3, 3)
    }

    alpha_up = (0.5 + theta_min) * gamma
    beta_up = (0.5 - theta_min) * gamma
    alpha_down = beta_up
    beta_down = alpha_up
    alpha_neu = (0.5 + theta_neu) * gamma
    beta_neu = (0.5 - theta_neu) * gamma
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha=c(alpha_neu,alpha_up,alpha_down), beta=c(beta_neu,beta_up,beta_down)),
        pn = list(size = DP),
        discrete = TRUE)

    hmm$states = states

    class(hmm) = 'dthmm.inhom'

    return(hmm)
}


#' Allele treeHMM with one theta state
#' @param pAD
#' @param DP
#' @param p_s
#' @param Q node marginals
#' @param Q_pair edge marginals
#' @param pa parent chain
#' @param ch children chains
#' @keywords internal
get_allele_treehmm = function(pAD, DP, p_s, Q = NULL, Q_pair = NULL, pa = NULL, ch = NULL, t = 1e-5, theta_min = 0.08, gamma = 20, prior = NULL) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('More than one gamma parameter')
    }
    
    # states
    states = c("neu", "theta_up", "theta_down")

    # no parent or children
    if (is.null(pa) & is.null(ch)) {
        calc_trans_mat = function(p_s, t) {
            matrix(
                c(
                    (1-t), t/2, t/2,
                    t, (1-t)*(1-p_s), (1-t)*p_s, 
                    t, (1-t)*p_s, (1-t)*(1-p_s)
                ),
                ncol = 3,
                byrow = TRUE
            )
        }
        As = lapply(p_s, function(p_s) {calc_trans_mat(p_s, t)})
        logAs = As %>% purrr::map(log)
    } else {
        logAs = treehmm_trans_mat(t, p_s, Q, Q_pair, pa, ch)
        As = logAs %>% purrr::map(exp)
    }
    
    # intitial probabilities
    if (is.null(prior)) {
        prior = rep(1/3, 3)
    }

    alpha_up = (0.5 + theta_min) * gamma
    beta_up = (0.5 - theta_min) * gamma
    alpha_down = beta_up
    beta_down = alpha_up
    alpha_neu = gamma/2
    beta_neu = gamma/2
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha=c(alpha_neu, alpha_up, alpha_down), beta=c(beta_neu, beta_up,beta_down)),
        pn = list(size = DP),
        discrete = TRUE)

    hmm$states = states
    hmm$logPi = logAs

    class(hmm) = 'dthmm.inhom'

    return(hmm)
}

# generate lookup table for conditional state probablities
# p_z[t, z_t_pa, z_t-1, z_t]
get_p_z = function(t, p_s) {

    # calculate state conditionals
    get_z_conditional = function(t, p_s, z_pa, z_t_1, z_t) {

        case_when(
            z_pa == 1 & z_t_1 == 1 & z_t == 1 ~ (1-t), 
            z_pa == 1 & z_t_1 == 1 & z_t == 2 ~ t/2, 
            z_pa == 1 & z_t_1 == 1 & z_t == 3 ~ t/2,
            z_pa == 1 & z_t_1 == 2 & z_t == 1 ~ t,
            z_pa == 1 & z_t_1 == 2 & z_t == 2 ~ (1-t)*(1-p_s),
            z_pa == 1 & z_t_1 == 2 & z_t == 3 ~ (1-t)*p_s,
            z_pa == 1 & z_t_1 == 3 & z_t == 1 ~ t, 
            z_pa == 1 & z_t_1 == 3 & z_t == 2 ~ (1-t)*p_s, 
            z_pa == 1 & z_t_1 == 3 & z_t == 3 ~ (1-t)*(1-p_s),
            # z_pa != 1 & z_t_1 == 1 & z_t == 1 ~ (1-t*eta),
            # z_pa != 1 & z_t_1 == 1 & z_t == 2 ~ t*eta/2, 
            # z_pa != 1 & z_t_1 == 1 & z_t == 3 ~ t*eta/2, 
            # z_pa != 1 & z_t_1 == 1 & z_t == 1 ~ t,
            # z_pa != 1 & z_t_1 == 1 & z_t == 2 ~ (1-t)/2, 
            # z_pa != 1 & z_t_1 == 1 & z_t == 3 ~ (1-t)/2, 
            z_pa != 1 & z_t_1 == 1 & z_t == 1 ~ t/2,
            z_pa != 1 & z_t_1 == 1 & z_t == 2 ~ t/2, 
            z_pa != 1 & z_t_1 == 1 & z_t == 3 ~ (1-t), 
            z_pa != 1 & z_t_1 == 2 & z_t == 1 ~ t,
            z_pa != 1 & z_t_1 == 2 & z_t == 2 ~ (1-t)*(1-p_s),
            z_pa != 1 & z_t_1 == 2 & z_t == 3 ~ (1-t)*p_s,
            z_pa != 1 & z_t_1 == 3 & z_t == 1 ~ t, 
            z_pa != 1 & z_t_1 == 3 & z_t == 2 ~ (1-t)*p_s, 
            z_pa != 1 & z_t_1 == 3 & z_t == 3 ~ (1-t)*(1-p_s)
        )
    }

    states = c("neu", "theta_up", "theta_down")

    states_grid = expand.grid(1:3, 1:3, 1:3) %>%
        setNames(c('z_t_pa', 'z_t_1', 'z_t'))

    p_z = sapply(
            states_grid %>% split(1:nrow(.)),
            function(Z){
                get_z_conditional(
                    t,
                    p_s = p_s,
                    Z$z_t_pa,
                    Z$z_t_1,
                    Z$z_t)
            }
        ) %>%
        array(dim = c(length(p_s),3,3,3))

    return(p_z)
}

# calculate transition matrices based on new marginals
#' @param pa parent index
#' @param ch children indices
#' @param Q Q[i, 1:t, z_t]
#' @param Q_pair Q_pair[i, t, z_t-1, z_t]
treehmm_trans_mat = function(t, p_s, Q, Q_pair, pa, ch) {

    bigT = length(p_s)

    eta = 1e2

    p_z = get_p_z(t, p_s)

    # defined for t=2:T
    logf_it = function(pa, ch, t, z_t_1, z_t) {

        if (!is.null(pa)) {
            # p_z[t, , z_t_1, z_t] is z_pa x z_t; Q[pa, t, ] is z_pa x 1 => m_pa is 1 x z_t
            m_pa = colSums(log(p_z[t, , z_t_1, z_t]) * Q[pa, t, ])
        } else {
            m_pa = log(p_z[t, 1, z_t_1, z_t])
        }

        if (!is.null(ch)) {
            m_ch = sapply(
                ch,
                function(i) {
                    sapply(
                        1:3,
                        function(w) {
                            sapply(
                                1:3,
                                function(y) {
                                    log(p_z[t, z_t, y, w]) * Q_pair[i, t, y, w]
                                }
                            ) %>% rowSums
                        }
                    ) %>% rowSums
                }
            ) %>%
            rowSums
            
        } else {
            m_ch = 0
        }

        m_pa + m_ch
    }

    # generate a list of transition matrices across time
    logf_z = lapply(
        1:bigT,
        function(t) {
            if (t == 1) {
                return(matrix(rep(NA,9), ncol = 3, nrow = 3))
            } else {
                sapply(1:3, function(z_t_1) {
                    logf_it(
                        pa = pa,
                        ch = ch,
                        t = t,
                        z_t_1 = z_t_1,
                        z_t = 1:3
                    )
                })
            }
        }
    )

    logf_z_old = logf_z

    # renormalize the transition matrix
    log_gz = rep(0,3)

    for (t in bigT:2) {
        logf_z[[t]] = t(t(logf_z[[t]]) + log_gz)
        log_gz = sapply(1:3, function(i){logSumExp(logf_z[[t]][i,])})
        logf_z[[t]] = logf_z[[t]] - log_gz
    }

    return(logf_z)
}

run_treehmm = function(bulks, pa_dict, ch_dict, theta_min = 0.08, gamma = 20, t = 1e-5, max_iter = 10) {
    
    bulks = bulks %>% split(.$sample)
    
    I = length(bulks)
    N = nrow(bulks[[1]])
    
    k = 1
    
    Q = array(rep(NA, I*N*3), dim = c(I, N, 3))
    Q_pair = array(rep(NA, I*N*3*3), dim = c(I, N, 3, 3))
    S = array(rep(NA, max_iter*I*N), dim = c(max_iter, I, N))
    P = array(rep(NA, max_iter*I*N), dim = c(max_iter, I, N))

    logprob = array(rep(NA, I*N*3), dim = c(I, N, 3))
    logtheta = log(get_p_z(t, bulks[[1]]$p_s))
    F = c()
    
    while (k <= max_iter) {

        message(glue('iter {k}'))
        
        ### update variational distribution ###
        for (i in I:1) {

            if (k == 1) {
                pa = NULL
                ch = NULL
            } else {
                pa = pa_dict[[i]]
                ch = ch_dict[[i]]
            }
            
            treehmm = bulks[[i]] %>% {
                get_allele_treehmm(
                    .$pAD, .$DP, .$p_s, 
                    t = t,
                    theta_min = theta_min,
                    Q = Q, 
                    gamma = gamma, 
                    Q_pair = Q_pair, 
                    pa = pa,
                    ch = ch
                    # prior = c(1-2*t, t, t)
                )
            }
            
            fb = forward_back_allele_R(treehmm)

            Q[i,,] = fb$marginals
            Q_pair[i,,,] = fb$edge_marginals
            logprob[i,,] = fb$logprob
            
            S[k,i,] = viterbi_allele(treehmm)
            P[k,i,] = Q[i,,1]
            
        }


        ### calculate free energy ###
        Fk = 0

        for (i in I:1) {

            pa = pa_dict[[i]]

            # entropy term
            H = sapply(
                2:N,
                function(t) {
                    q_cond = Q_pair[i, t, , ]/Q[i, t-1, ]
                    sum(Q_pair[i, t, , ] * log(q_cond))
                }
            ) %>% sum

            H = H + sum(Q[i,1,] * log(Q[i,1,]))

            # likelihood term
            L = sapply(
                2:N,
                function(t) {
                    
                    emission_term = sum(Q[i,t,] * logprob[i,t,])
                    
                    transition_term = sum(sapply(1:3,
                        function(z_pa) {
                            sum(logtheta[t, z_pa, , ] * Q_pair[i, t, , ] * Q[pa, t, z_pa])
                        })
                    )
                    
                    emission_term + transition_term
                }
            ) %>% sum

            # free energy
            Fk = Fk + H - L

        }

        F = c(F, Fk)

        k = k + 1
    }

    return(list(S = S, P = P, F = F))
}



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
        
        l_x = p_x(x = x,  getj(obj$pm, k), obj$pn, log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })
        
    logphi <- log(as.double(obj$delta))

    if (is.null(obj$logPi)) {
        logPi <- lapply(obj$Pi, log)
    } else {
        logPi = obj$logPi
    }
            
    marginals = forward_backward_compute(obj, logphi, logprob, logPi, n, m)

    colnames(marginals) = obj$states

    return(marginals)
}

forward_back_allele_R = function (obj, ...) {

    # case of one-data point
    if (length(obj$x) == 1) {
        return(NA)
    }
    
    x <- obj$x
    p_x <- makedensity(obj$distn)
    
    m <- nrow(obj$Pi[[1]])
    n <- length(x)
    
    logprob = sapply(1:m, function(k) {
        
        l_x = p_x(x = x,  getj(obj$pm, k), obj$pn, log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })
        
    logphi <- log(as.double(obj$delta))
    logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
    lscale <- as.double(0)

    if (is.null(obj$logPi)) {
        logPi <- lapply(obj$Pi, log)
    } else {
        logPi = obj$logPi
    }

    for (t in 1:n) {
        
        if (t > 1) {
            logphi <- sapply(1:m, function(j) matrixStats::logSumExp(logphi + logPi[[t]][,j]))
        }
                          
        logphi <- logphi + logprob[t,]
                          
        logsumphi <- matrixStats::logSumExp(logphi)
                          
        logphi <- logphi - logsumphi
                          
        lscale <- lscale + logsumphi
                          
        logalpha[t,] <- logphi + lscale
                          
        LL <- lscale
    }

    logbeta <- matrix(as.double(rep(0, m * n)), nrow = n)
    logphi <- log(as.double(rep(1/m, m)))
    lscale <- as.double(log(m))

    for (t in seq(n-1, 1, -1)){
        
        logphi = sapply(1:m, function(i) matrixStats::logSumExp(logphi + logprob[t+1,] + logPi[[t+1]][i,]))

        logbeta[t,] <- logphi + lscale

        logsumphi <- matrixStats::logSumExp(logphi)

        logphi <- logphi - logsumphi

        lscale <- lscale + logsumphi
    }

    marginals = exp(logalpha + logbeta - LL)
    colnames(marginals) = obj$states

    edge_marginals = lapply(
        1:(n-1),
        function(t) {
            x = exp(outer(logalpha[t,], logbeta[t+1,], FUN = '+') + logPi[[t+1]] + logprob[t+1,] - LL)
            return(x)
        }
    )
    # defined for (t-1, t) for t=2:T
    edge_marginals = c(matrix(rep(NA,9), ncol = 3, nrow = 3), edge_marginals)

    edge_marginals = aperm(array(unlist(edge_marginals), dim = c(m,m,n), dimnames = list(obj$states, obj$states, 1:(n))), c(3,1,2))

    return(list(marginals = marginals, edge_marginals = edge_marginals, logalpha = logalpha, logbeta = logbeta, logprob = logprob, LL = LL, logPi = logPi))
}


# only compute total log likelihood
#' @keywords internal
likelihood_allele = function (obj, ...) {
        
    x <- obj$x
    p_x <- makedensity(obj$distn)
    
    m <- nrow(obj$Pi[[1]])
    n <- length(x)
    
    logprob = sapply(1:m, function(k) {
        
        l_x = p_x(x = x,  getj(obj$pm, k), obj$pn, log = TRUE)
        
        l_x[is.na(l_x)] = 0
        
        return(l_x)
        
    })

    logphi <- log(as.double(obj$delta))
    
    logPi <- lapply(obj$Pi, log)
        
    LL <- likelihood_allele_compute(obj, logphi, logprob, logPi, n, m)

    return(LL)
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
            beta = c(beta_up, beta_down)), pn = list(size = DP), 
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

#' @keywords internal
# calc_allele_maxlik = function (pAD, DP, p_s, theta, gamma = 20) {
#     hmm = get_allele_hmm(pAD, DP, p_s, theta, gamma)
#     LL = viterbi_allele(hmm)$LL
#     return(LL)
# }

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
        nu[1, ] = nu[1, ] + dfunc(x=x[1], object$pm, getj(object$pn, 1), log = TRUE)
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
            nu[i, ] = nu[i, ] + dfunc(x=x[i], object$pm, getj(object$pn, i), log = TRUE)
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
        
    return(list(y = y, LL = LL))
}

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

#' @keywords internal
calc_trans_mat = function(t, p_s, w, states_cn, states_phase) {

    sapply(1:length(states_cn), function(from) {
        sapply(1:length(states_cn), function(to) {
            get_trans_probs(t, p_s, w, states_cn[from], states_phase[from], states_cn[to], states_phase[to])
        }) %>% t
    }) %>% t %>%
    array(dim = c(length(states_cn), length(states_cn), length(p_s)))

}

########## HMM wrappers ###########
#' @keywords internal
run_hmm_mv_inhom = function(
    pAD, DP, p_s, Y_obs = 0, lambda_ref = 0, d_total = 0, theta_min = 0.08, theta_neu = 0,
    bal_cnv = TRUE, phi_del = 2^(-0.25), phi_amp = 2^(0.25), phi_bamp = phi_amp, phi_bdel = phi_del, 
    alpha = 1, beta = 1, 
    mu = 0, sig = 1,
    exp_model = 'gpois',
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
        pn = list(size = DP),
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




