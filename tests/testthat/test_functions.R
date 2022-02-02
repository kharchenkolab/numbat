
library(Numbat)

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



