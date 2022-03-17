
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
      gtf = gtf_hg38,
      genetic_map = genetic_map_hg38
  )

  expect_true(length(bulk) > 0)

  message('running HMM..')

  bulk = analyze_bulk(bulk, t = 1e-5)

  expect_true(length(bulk) > 0)
  
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

  LL = likelihood_allele(pre_likelihood_hmm)
  expect_equal(as.integer(LL), -736)

})



test_that("Check that forward_backward() works as expected", {

  p_major = forward_back_allele(pre_likelihood_hmm)
  expect_equal(is.vector(p_major), TRUE)
  expect_equal(length(p_major), 1042)
  expect_equal(round(p_major[1], 3), 0.963)
  expect_equal(round(p_major[2], 3), 0)
  expect_equal(round(p_major[3], 3), 0)
  expect_equal(round(p_major[10], 3), 0.745)

})

