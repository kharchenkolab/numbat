
library(numbat)

test_that("dpoilog works", {
  expect_equal(dpoilog(c(1,11),c(1,1),c(1,1)), c(0.175733342664327, 0.0150105250670325))
})

test_that("HMM works", {

  message('getting pseudobulk ..')

  bulk = get_bulk(
      count_mat = count_mat_ATC2,
      df_allele = df_allele_ATC2,
      lambdas_ref = ref_hca,
      gtf_transcript = gtf_hg38,
      genetic_map = genetic_map_hg38
  )

  expect_equal(length(bulk), 28)

  message('running HMM..')

  bulk = analyze_bulk(bulk, t = 1e-5)

  expect_equal(length(bulk), 78)
  
})



test_that("logSumExp() works", {

  a = matrixStats::logSumExp(c(1.2, 3.4))
  b = numbat:::logSumExp(c(1.2, 3.4))
  expect_equal(a, b)
  expect_equal(a, 3.5050833)
  expect_equal(b, 3.5050833)  
  c = matrixStats::logSumExp(c(1.2, 3.4, -5.6, -7.8))
  d = numbat:::logSumExp(c(1.2, 3.4, -5.6, -7.8))
  expect_equal(c, d)  
  expect_equal(c, 3.5052067)
  expect_equal(d, 3.5052067) 
   
})

test_that("Check that likelihood_allele() works as expected", {

  x <- pre_likelihood_hmm$x
  p_x <- numbat:::makedensity(pre_likelihood_hmm$distn)
  m <- nrow(pre_likelihood_hmm$Pi[[1]])
  n <- length(x)
  logprob = sapply(1:m, function(k) {
      l_x = p_x(x = x,  numbat:::getj(pre_likelihood_hmm$pm, k), pre_likelihood_hmm$pn, log = TRUE)
      l_x[is.na(l_x)] = 0
      return(l_x)
  })
  logphi <- log(as.double(pre_likelihood_hmm$delta))
  logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
  LL <- as.double(0)
  logPi <- lapply(pre_likelihood_hmm$Pi, log)

  LL <- numbat:::likelihood_allele_compute(pre_likelihood_hmm, logphi, logprob, logPi, n, m)
  expect_equal(as.integer(LL), -736)

})


