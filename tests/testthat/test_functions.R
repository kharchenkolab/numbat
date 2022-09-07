
library(numbat)

test_that("dpoilog works", {
  expect_equal(numbat:::dpoilog(c(1,11),c(1,1),c(1,1)), c(0.175733342664327, 0.0150105250670325))
})

# test_that("HMM works", {

#   message('getting pseudobulk ..')

#   bulk = get_bulk(
#       count_mat = count_mat_ATC2,
#       df_allele = df_allele_ATC2,
#       lambdas_ref = ref_hca,
#       gtf = gtf_hg38
#   )

#   expect_true(length(bulk) > 0)

#   message('running HMM..')

#   bulk = analyze_bulk(bulk, t = 1e-5)

#   expect_true(length(bulk) > 0)
  
# })

test_that("logSumExp() works", {

  b = numbat:::logSumExp(c(1.2, 3.4))
  expect_equal(b, 3.5050833)  
  d = numbat:::logSumExp(c(1.2, 3.4, -5.6, -7.8))
  expect_equal(d, 3.5052067) 
   
})

test_that("Check that likelihood_allele() works as expected", {

  LL = numbat:::likelihood_allele(pre_likelihood_hmm)
  expect_equal(as.integer(LL), -736)

})

test_that("Check that forward_backward() works as expected", {

  p_major = numbat:::forward_back_allele(pre_likelihood_hmm)[,1]
  expect_equal(is.vector(p_major), TRUE)
  expect_equal(length(p_major), 1042)
  expect_equal(round(p_major[1], 3), 0.963)
  expect_equal(round(p_major[2], 3), 0)
  expect_equal(round(p_major[3], 3), 0)
  expect_equal(round(p_major[10], 3), 0.745)

})

test_that("Check that viterbi_allele() works as expected", {

  states = numbat:::viterbi_allele(pre_likelihood_hmm)
  expect_equal(sum(states == 1), 440)
  expect_equal(sum(states == 2), 602)
  expect_equal(states[1], 1)
  expect_equal(states[2], 2)
  expect_equal(states[3], 2)
  expect_equal(states[1024], 1)
  
})



test_that("Check that roman2int works as expected", {
   val.3899 <- 'MMMDCCCXCIX'
   val.3900 <- 'MMMCM'
   val.4000 <- 'MMMM'
   expect_equal(names(roman2int(val.3899)), "MMMDCCCXCIX")
   expect_equal(as.integer(roman2int(val.3899)), 3899)
   expect_equal(names(roman2int(val.3900)), "MMMCM")
   expect_equal(as.integer(roman2int(val.3900)), 3900)
   expect_equal(names(roman2int(val.4000)), "MMMM")
   expect_equal(as.integer(roman2int(val.4000)), 4000)   
})
