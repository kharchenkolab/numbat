
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


test_that("Check that likelihood_allele() works as expected", {

  x <- pre_likelihood_hmm$x
  p_x <- makedensity(pre_likelihood_hmm$distn)
  m <- nrow(pre_likelihood_hmm$Pi[[1]])
  n <- length(x)
  logprob = sapply(1:m, function(k) {
      l_x = p_x(x = x,  getj(pre_likelihood_hmm$pm, k), pre_likelihood_hmm$pn, log = TRUE)
      l_x[is.na(l_x)] = 0
      return(l_x)
  })
  logphi <- log(as.double(pre_likelihood_hmm$delta))
  logalpha <- matrix(as.double(rep(0, m * n)), nrow = n)
  LL <- as.double(0)
  logPi <- lapply(pre_likelihood_hmm$Pi, log)

  LL <- numbat:::likelihood_allele_compute(obj, logphi, logprob, logPi, n, m)
  expect_equal(as.integer(LL) -736)

})


