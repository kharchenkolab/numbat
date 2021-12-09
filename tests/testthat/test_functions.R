test_that("dpoilog works", {
  expect_equal(dpoilog(c(1,11),c(1,1),c(1,1)), c(0.175733342664327, 0.0150105250670325))
})

test_that("HMM works", {

  message('getting pseudobulk ..')

  bulk = get_bulk(
      count_mat = count_mat_ATC2,
      df = df_allele_ATC2,
      lambdas_ref = ref_hca,
      gtf_transcript = gtf_transcript,
      genetic_map = genetic_map_hg38,
      min_depth = 0
  )

  message('running HMM..')

  bulk = bulk %>% analyze_bulk_lnpois(t = 1e-5)
  
})


